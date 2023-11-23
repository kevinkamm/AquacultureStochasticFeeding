function [params,aFin,resOuter,varargout]=cortazarDaily1C(dataName,dataStartDate,dataEndDate,r,varargin)

verbose=0;
for k=1:length(varargin)
    switch varargin{k}
        case 'verbose'
            verbose = varargin{k+1};
    end
end

%% Extract data from spreadsheet
disp('Load data')
% dataName='Soy';
data=load(['Data/Daily/',dataName]);
prices=data.output;
dates=data.udates;
ttm=data.uTTM';

filledPrices=fillmissing(prices,'nearest',2);

% dataStartDate=datetime(2018,1,1);
% dataEndDate=datetime(2022,12,31);
dateInd = dates >=dataStartDate & dates <=dataEndDate;

yData=log(prices(dateInd,:));
% yData=log(filledPrices(dateInd,:));

%%
disp('Calibrate')
ticCal=tic;
[params,resOuter,aFin,resInner]=calibrate(r,ttm,yData);
ctimeCal=toc(ticCal);
fprintf('Calibration done in %g s with outer error %g and inner mean error %g\n',ctimeCal,resOuter,mean(resInner));

%%
if verbose>0
    [sigma11,sigma12,kappa1,alpha1,lambda1,rho,paramDesc]=paramUnpack(params);
    paramTable=array2table(reshape(params,1,[]),'VariableNames',paramDesc);
    disp(paramTable)
    disp(join(split(num2str([nan,params(1:6),aFin(end,2),aFin(end,1)]),whitespacePattern),', '))
end
%%
if verbose > 1
    fig=newFigure();
    plot(filledPrices(dateInd,1),'bx')
    plot(aFin(:,1),'r--')
end
if nargout > 3
    varargout{1}=dates(dateInd);
end

end
function [sigma11,sigma12,kappa1,alpha1,lambda1,rho,varargout]=paramUnpack(params)
    sigma11=params(1);
    sigma12=params(2);
    kappa1=params(3);
    alpha1=params(4);
    lambda1=params(5);
    rho=params(6);
    if nargout > 6
        paramDesc={'sigma1','sigma2','kappa','alpha','lambda','rho'};
        varargout{1}=paramDesc;
    end
end

function [pFin,resOuter,aFin,resInner]=calibrate(r,ttm,yData)

%
% dimensions: weeks x time-to-maturity
%

FMarket=yData;% log futures
indNan=~isnan(yData);
[nDates,nTTM] = size(yData);
FMarket(~indNan)=0;

optionsOuter=optimoptions('fmincon','Display','final-detailed');
% optionsOuter=optimoptions('lsqnonlin','MaxFunctionEvaluations',1e4,'Display','final-detailed');
optionsInner=optimset('Display','notify');

p0=[0.01, 0.01, 0.01, 0.01, 0.01, 0.9];
% p0=[0.9, 0.9, 0.9, 0.9, 0.9, 0.99];
lbP=1e-4.*ones(6,1);
lbP(end)=-1;
lbP(4)=-4;
ubP=4.*ones(6,1);ubP(end)=1;

% a01=yData(1,:);
% a01=a01(~isnan(a01));
% a01=a01(1);
% a0=[a01,r];

[pFin,resOuter,aFin,resInner]=outerCalibration(p0);
aFin(:,1)=exp(aFin(:,1));

function [p,resOuter,aCal,errInner]=outerCalibration(p0)
    % p0=ga(@(p)mean(abs(outerObjective(a0,p)).^2,'all'),length(p0),[],[],[],[],lbP,ubP,[],[],optionsOuterGA)
    % [p,resOuter]=lsqnonlin(@(p)outerObjective(p),p0,lbP,ubP,optionsOuter);
    [p,resOuter]=fmincon(@(p)outerObjective(p),p0,[],[],[],[],lbP,ubP,[],optionsOuter);

    [aCal,errInner]=innerCalibration(p);
end

function res=outerObjective(p)
    [aCal,innerErr]=innerCalibration(p);
    Fmodel=logFutures(aCal,p);
    res=Fmodel(indNan)-FMarket(indNan);
    % res=mean(abs(res).^2,'all');
    res=mean(abs(res).^2,'all')+mean(innerErr,'all');
end

function [aCal,errInner]=innerCalibration(p)
    aCal=zeros(nDates,2);
    errInner=zeros(nDates,1);
    for ti=1:1:nDates
        % [a,res]=fminunc(@(a)innerObjective(a,p,ti),a0,optionsInner);
        % a0=a;
        % aCal(ti,:)=exp(a0);

        [A,b]=logFuturesLeastSquare(p,ti);
        % [a,res]=lsqnonneg(A,b,optionsInner);

        a=A\b; res=norm(A*a-b);

        aCal(ti,1)=a(1);
        aCal(ti,2)=a(2);

        errInner(ti)=res;
    end
    
end

% function res=innerObjective(a,p,ti)
%     Fmodel=logFutures(exp(a),p);
%     % res=mean(abs(Fmodel(1)-FMarket(ti,1)).^2,'all');%changed here for spot to smallest maturity not all of them
%     res=mean(abs(Fmodel-FMarket(ti,:)).^2,'all');%changed here for spot to smallest maturity not all of them
% end
% 
function Fmodel=logFutures(a,p)

    [sigma11,sigma12,kappa1,alpha1,lambda1,rho12]=paramUnpack(p);
    
    tau1=ttm;

    Atau1=(r-alpha1+lambda1./kappa1 + sigma12.^2./(2.*kappa1.^2)-sigma11.*sigma12.*rho12./kappa1).*tau1+...
    sigma12.^2.*((1-exp(-2.*kappa1.*tau1))./(kappa1.^3) )./4+...
    (alpha1.*kappa1-lambda1+sigma11.*sigma12.*rho12-sigma12.^2./kappa1).*((1-exp(-kappa1.*tau1))./(kappa1.^2));

    Fmodel=a(:,1)-a(:,2).*(1-exp(-kappa1.*tau1))./kappa1 + Atau1;
end

function [A,b]=logFuturesLeastSquare(p,ti)
    % A, b for:
    %   min_{x>0} \norm{A x - b} (lsqnonneg)
    % x\in R^d, b\in R^m: A\in R^{m,d}
    % Here x\in R^2, b=Fmarket - Atau

    [sigma11,sigma12,kappa1,alpha1,lambda1,rho12]=paramUnpack(p);
    
    % indTTM = ~isnan(yData(ti,:));
    FM=yData(ti,indNan(ti,:));
    tau1=ttm(indNan(ti,:));

    Atau1=(r-alpha1+lambda1./kappa1 + sigma12.^2./(2.*kappa1.^2)-sigma11.*sigma12.*rho12./kappa1).*tau1+...
    sigma12.^2.*((1-exp(-2.*kappa1.*tau1))./(kappa1.^3) )./4+...
    (alpha1.*kappa1-lambda1+sigma11.*sigma12.*rho12-sigma12.^2./kappa1).*((1-exp(-kappa1.*tau1))./(kappa1.^2));
    
    A=[ones(length(tau1),1),-(1-exp(-kappa1.*tau1'))./kappa1];
    b=(FM-Atau1)';

end

end
