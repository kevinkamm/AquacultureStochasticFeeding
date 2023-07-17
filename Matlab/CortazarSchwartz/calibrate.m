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