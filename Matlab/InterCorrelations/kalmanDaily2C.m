function [params,ss_att,negLogLikelihood]=kalmanDaily2C(dataStartDate,dataEndDate,dataPanel,r,varargin)

verbose=0;
for k=1:length(varargin)
    switch varargin{k}
        case 'verbose'
            verbose = varargin{k+1};
    end
end

% dataStartDate=datetime(2018,1,1);
% dataEndDate=datetime(2022,12,31);
% % dataPanel=[[0,1,2,3,4,5]',[1,2,3,4,5,6]']./12; % small panel
% % dataPanel=[[0,1,5,6,10,11]',1+[0,1,5,6,10,11]']./12; % mixed panel
% % dataPanel=[[5,6,10,11,14,15]',1+[5,6,10,11,14,15]']./12; % medium panel
% dataPanel=[[10,11,14,15,19,20]',1+[10,11,14,15,19,20]']./12; % large panel


%% Salmon Data
verb(@()disp('Load data'),verbose)
panelSalmon=zeros(1,size(dataPanel,1));

dataName='Salmon';
data=load(['Data/Daily/',dataName]);
prices=data.output;
dates=data.udates;
ttm=data.uTTM;


dateInd = dates >=dataStartDate & dates <=dataEndDate;
ttmInd = ttm >= 0 & ttm <= dataPanel(end);

filledPrices=fillmissing(prices,'nearest',2);
ySalmon=zeros(sum(dateInd),size(dataPanel,1));
for p=1:size(dataPanel,1)
    indP= ttm>dataPanel(p,1) & ttm<=dataPanel(p,2);
    ySalmon(:,p)=mean(filledPrices(dateInd,indP),2);
    panelSalmon(1,p)=mean(ttm(indP));
end

datesSalmon=dates(dateInd);
ySalmon=log(ySalmon);
indContracts=1:size(ySalmon,2); % start from 2 to exclude spot
panelSalmon=panelSalmon(:,indContracts);

%% Soy Data
verb(@()disp('Load data'),verbose)
panelSoy=zeros(1,size(dataPanel,1));

dataName='Soy';
data=load(['Data/Daily/',dataName]);
prices=data.output;
dates=data.udates;
ttm=data.uTTM;


dateInd = dates >=dataStartDate & dates <=dataEndDate;
ttmInd = ttm >= 0 & ttm <= dataPanel(end);

filledPrices=fillmissing(prices,'nearest',2);
ySoy=zeros(sum(dateInd),size(dataPanel,1));
for p=1:size(dataPanel,1)
    indP= ttm>dataPanel(p,1) & ttm<=dataPanel(p,2);
    ySoy(:,p)=mean(filledPrices(dateInd,indP),2);
    panelSoy(1,p)=mean(ttm(indP));
end

datesSoy=dates(dateInd);
ySoy=log(ySoy);
indContracts=1:size(ySoy,2); % start from 2 to exclude spot
panelSoy=panelSoy(:,indContracts);


%% Joint Data
[C,iSoy,iSalmon]=intersect(datesSoy',datesSalmon');

% dt = 1/365; % daily time step
dt = years(mean(diff(C))); % =250 daily time step

ySalmon=ySalmon(iSalmon,:);
ySoy=ySoy(iSoy,:);

y=[ySalmon,ySoy];

% panel=[panelSalmon(indContractsSalmon),panelSoy(indContractsSoy)];
panel1=panelSalmon;
panel2=panelSoy;

a0=[mean(ySalmon(:,1));r;mean(ySoy(:,1));r];
P0=zeros(4,4);

%%
disp('Calibrate')
verb(@()disp('Calibrate'),verbose)
% varEpsMode='single';
varEpsMode='diag';
% varEpsMode='tridiag';
ticEM=tic;
[params,negLogLikelihood]=EMschwartz(y,r,dt,panel1,panel2,a0,P0,varEpsMode);
ctimeEM=toc(ticEM);
verb(@()fprintf('Elapsed time for EM-Algo: %g s with log-likelihood %g\n',ctimeEM,-negLogLikelihood),verbose)

[P0,a0,c,T,varEta,dtau,Ztau,varEpsMode,~]=schwartzSSM(params,y,r,dt,panel1,panel2,a0,P0,varEpsMode);
[logL,ss_att]=kalmanFilter(P0,a0,c,T,varEta,dtau,Ztau,varEpsMode,y);



%%
if verbose > 0
    %%
    [mu1,sigma11,sigma12,kappa1,alpha1,lambda1,mu2,sigma21,sigma22,kappa2,alpha2,lambda2,rho12,rho13,rho14,rho23,rho24,rho34,s,paramDesc]=paramUnpack(params);
    paramTable=array2table(params,'VariableNames',paramDesc);
    disp(paramTable)
    
    %%
    if verbose >1
        fig=newFigure();
        subplot(1,2,1);hold on;
        plot(exp(ySalmon(:,1)),'k','linewidth',1);%spot
        plot(exp(ySalmon(:,2:end)),'r--','linewidth',1);%futures
        plot(exp(ss_att(:,1)),'b','linewidth',1);%Kalman
        title('Salmon')
        
        subplot(1,2,2);hold on;
        plot(exp(ySoy(:,1)),'k','linewidth',1);%spot
        plot(exp(ySoy(:,2:end)),'r--','linewidth',1);%futures
        plot(exp(ss_att(:,3)),'b','linewidth',1);%Kalman
        title('Soy')
    end
    %%
    disp('Salmon: mu1, sigma11, sigma12, kappa1, alpha1, lambda1, delta10, Price10')
    disp(join(split(num2str([params(1:6),ss_att(1,2),exp(ss_att(1,1))]),whitespacePattern),', '))
    disp('Soy: mu2, sigma21, sigma22, kappa2, alpha2, lambda2, delta20, Price20')
    disp(join(split(num2str([params(7:12),ss_att(1,4),exp(ss_att(1,3))]),whitespacePattern),', '))
    disp('Correlations: rho12, rho13, rho14, rho23, rho24, rho34')
    disp(join(split(num2str(params(13:18)),whitespacePattern),', '))
end


end
function verb(func,flag)
    if flag>0
            func();
    end
end

function [params,negLogLikelihood]=EMschwartz(y,r,dt,panel1,panel2,a0,P0,varEpsMode)
    [nObservations,nContracts] = size(y);
    
    switch varEpsMode
        case 'single'
            n=1;
        case 'diag'
            n=nContracts;
        case 'tridiag'
            n=nContracts*(nContracts+1)/2;
        otherwise
            error('Unrecognized mode for variance of epsilon')
    end
    
    %params=[mu1,sigma11,sigma12,kappa1,alpha1,lambda1,...
    %        mu2,sigma21,sigma22,kappa2,alpha2,lambda2,...
    %        rho12,rho13,rho14,rho23,rho24,rho34,...
    %        s1,...,sn];
    %           

    lb=1e-2.*ones(1,6*2+6+n);
    lb(12:18)=-1; %correlations
    ub=5.*ones(1,6*2+6+n);
    ub(12:end)=1; %correlations & s
    lb(5)=-ub(5);%alpha1
    lb(11)=-ub(11);%alpha2

    x0=.1.*ones(1,6*2+6+n);
    x0(12:18)=.5; %correlations

    optionsF = optimoptions('fmincon','Display','final-detailed','UseParallel',false,'StepTolerance',1e-12,'OptimalityTolerance',1e-9,'MaxFunctionEvaluations',1e6,'MaxIterations',1e6);

    [params,negLogLikelihood]= fmincon(...
        @(params) objective(params,y,r,dt,panel1,panel2,a0,P0,varEpsMode),...
        x0,[],[],[],[],lb,ub,[],optionsF);

end
function logL=objective(params,y,r,dt,panel1,panel2,a0,P0,varEpsMode)
    [P0,a0,c,T,varEta,dtau,Ztau,varEps,y]=schwartzSSM(params,y,r,dt,panel1,panel2,a0,P0,varEpsMode);
    logL=kalmanFilter(P0,a0,c,T,varEta,dtau,Ztau,varEps,y);
end


function [P0,a0,c,T,varEta,dtau,Ztau,varEps,y]=schwartzSSM(params,y,r,dt,panel1,panel2,a0,P0,varEpsMode)
%%SCHWARTZSSM is a parameter map for the Schwartz-two-factor SSM for the
% Matlab built-in SSM.
%
%   Schwartz-2-factor model:
%   dP_t = (r-\delta_t) P_t dt + \sigma_1 P_t dW^1_t
%   d\delta_t = (\kappa (\alpha-\delta_t)-\lambda) dt + \sigma_2 dW^2_t
%   d<W^1,W^2>_t = \rho dt
%
% tau is the maturities in the panel (could be a vector).
%
% Model SSM:
%   State equation: [x(t+1);delta(t+1)] = c + T [x(t);delta(t)] + eta
%   Observation equation: y(t+1,tau) = d(tau) + Z(tau) [x(t);delta(t)]  + epsilon
%   Var(eta) = [sigma_1^2,              rho sigma_1 sigma_2; 
%               rho sigma_1 sigma_2,               sigma_2^2]*dt
%   c=[(mu-sigma_1^2/2);kappa alpha] * dt
%   T = [1, -dt; 0, 1- kappa dt]
%   d(tau) = A(tau)
%   Z(tau) = [1,-(1-e^(-kappa tau))/kappa]
%
%
    panel1=panel1(:); % time to maturities for different contracts in panel1
    panel2=panel2(:); % time to maturities for different contracts in panel
    [mu1,sigma11,sigma12,kappa1,alpha1,lambda1,mu2,sigma21,sigma22,kappa2,alpha2,lambda2,rho12,rho13,rho14,rho23,rho24,rho34,s]=paramUnpack(params);

    d=length(panel1)+length(panel2);
    
    switch varEpsMode
        case 'single'
            varEps= s.^2.*eye(d);
        case 'diag'
            varEps= diag(s.^2);
        case 'tridiag'
            tmp=reshape(1:d^2,d,d);
            ind=triu(tmp);ind=ind(ind>0);
            tmp=zeros(d,d);
            tmp(ind)=s;
            varEps=tmp'*tmp;
    end


    % [nObservations,nContracts] = size(y);
    
    % cov11=-(2.*rho12.*sigma12.*sigma11.*exp(kappa1.*(-dt)).*(exp(kappa1.*dt).*(kappa1.*dt-1)+1))./(kappa1.^2)+(sigma12.^2.*exp(-2.*kappa1.*dt).*(exp(2.*kappa1.*dt).*(2.*kappa1.*dt-3)+4.*exp(kappa1.*dt)-1))./(2.*kappa1.^3)+sigma11.^2.*dt;
    % cov12=(sigma12.*exp(-2.*kappa1.*dt).*(exp(kappa1.*dt)-1).*(sigma12+exp(kappa1.*dt).*(2.*kappa1.*rho12.*sigma11-sigma12)))./(2.*kappa1.^2);
    % cov13=(exp(-2.*(kappa1+kappa2).*dt).*(kappa1.^3.*sigma11.*exp((2.*kappa1+kappa2).*dt).*(exp(kappa2.*dt).*(kappa2.^2.*rho13.*sigma21.*dt+rho14.*sigma22.*(1-kappa2.*dt))-rho14.*sigma22)+kappa1.^2.*exp((2.*kappa1+kappa2).*dt).*(sigma22.*(rho24.*sigma12-kappa2.*rho14.*sigma11)+exp(kappa2.*dt).*(kappa2.^2.*sigma21.*dt.*(kappa2.*rho13.*sigma11-rho23.*sigma12)-sigma22.*(kappa2.*dt-1).*(kappa2.*rho14.*sigma11-rho24.*sigma12)))+kappa2.*kappa1.*sigma12.*exp((kappa1+kappa2).*dt).*(kappa2.^2.*rho23.*sigma21.*dt.*(-exp((kappa1+kappa2).*dt))+kappa2.*exp(kappa2.*dt).*(rho23.*sigma21.*(exp(kappa1.*dt)-1)+rho24.*sigma22.*dt.*exp(kappa1.*dt))-rho24.*sigma22.*(exp(kappa1.*dt)-1).*(exp(kappa2.*dt)-1))+kappa2.^2.*sigma12.*exp((kappa1+2.*kappa2).*dt).*(exp(kappa1.*dt)-1).*(kappa2.*rho23.*sigma21-rho24.*sigma22)))./(kappa1.^2.*kappa2.^2.*(kappa1+kappa2));
    % cov14=(sigma22.*exp((kappa1+kappa2).*(-dt)).*(kappa1.^2.*rho14.*sigma11.*exp(kappa1.*dt).*(exp(kappa2.*dt)-1)+kappa1.*exp(kappa1.*dt).*(exp(kappa2.*dt)-1).*(kappa2.*rho14.*sigma11-rho24.*sigma12)+kappa2.*rho24.*sigma12.*(exp(kappa1.*dt)-1)))./(kappa1.*kappa2.*(kappa1+kappa2));
    % cov22=(sigma12.^2.*(1-exp(-2.*kappa1.*dt)))./(2.*kappa1);
    % cov23=(sigma12.*exp((kappa1+kappa2).*(-dt)).*(kappa2.*exp(kappa2.*dt).*(exp(kappa1.*dt)-1).*(kappa2.*rho23.*sigma21-rho24.*sigma22)+kappa1.*(kappa2.*rho23.*sigma21.*exp(kappa2.*dt).*(exp(kappa1.*dt)-1)+rho24.*sigma22.*(exp(kappa2.*dt)-1))))./(kappa1.*kappa2.*(kappa1+kappa2));
    % cov24=(rho24.*sigma12.*sigma22.*exp((kappa1+kappa2).*(-dt)).*(exp((kappa1+kappa2).*dt)-1))./(kappa1+kappa2);
    % cov33=-(2.*rho34.*sigma22.*sigma21.*exp(kappa2.*(-dt)).*(exp(kappa2.*dt).*(kappa2.*dt-1)+1))./(kappa2.^2)+(sigma22.^2.*exp(-2.*kappa2.*dt).*(exp(2.*kappa2.*dt).*(2.*kappa2.*dt-3)+4.*exp(kappa2.*dt)-1))./(2.*kappa2.^3)+sigma21.^2.*dt;
    % cov34=(sigma22.*exp(-2.*kappa2.*dt).*(exp(kappa2.*dt)-1).*(sigma22+exp(kappa2.*dt).*(2.*kappa2.*rho34.*sigma21-sigma22)))./(2.*kappa2.^2);
    % cov44=(sigma22.^2.*(1-exp(-2.*kappa2.*dt)))./(2.*kappa2);

    cov11 = sigma11.^2;
    cov12 = sigma11 .* sigma12 .* rho12;
    cov13 = sigma11 .* sigma21 .* rho13;
    cov14 = sigma11 .* sigma22 .* rho14;
    cov22 = sigma12.^2;
    cov23 = sigma12 .* sigma21 .* rho23;
    cov24 = sigma12 .* sigma22 .* rho24;
    cov33 = sigma21.^2;
    cov34 = sigma21 .* sigma22 .* rho34;
    cov44 = sigma22.^2;

    varEta = [cov11,cov12,cov13,cov14;
              cov12,cov22,cov23,cov24;
              cov13,cov23,cov33,cov34;
              cov14,cov24,cov34,cov44].*dt;

    c=[mu1-sigma11./2;
        kappa1.*alpha1;
        mu2-sigma21./2;
        kappa2.*alpha2].*dt;

    T = [1  -dt             0               0;
         0  1-kappa1.*dt    0               0;
         0  0               1               -dt;
         0  0               0               1-kappa2.*dt];

    % c=[(mu1-sigma11^2/2-alpha1).*dt + (1-exp(-kappa1.*dt)).*alpha1./kappa1;
    %     alpha1.*(1-exp(-kappa1.*dt));
    %     (mu2-sigma21^2/2-alpha2).*dt + (1-exp(-kappa2.*dt)).*alpha2./kappa2;
    %     alpha2.*(1-exp(-kappa2.*dt))];
    % T=[1, (exp(-kappa1.*dt)-1)./kappa1, 0, 0; 
    %    0, exp(-kappa1.*dt), 0, 0;
    %    0, 0, 1, (exp(-kappa2.*dt)-1)./kappa2;
    %    0, 0, 0, exp(-kappa2.*dt)];
    
    dtau1 = (r-alpha1+lambda1./kappa1 + sigma12.^2./(2.*kappa1.^2)-sigma11.*sigma12.*rho12./kappa1).*panel1+...
    sigma12.^2.*((1-exp(-2.*kappa1.*panel1))./(kappa1.^3) )./4+...
    (alpha1.*kappa1-lambda1+sigma11.*sigma12.*rho12-sigma12.^2./kappa1).*((1-exp(-kappa1.*panel1))./(kappa1.^2));
    
    Ztau1=[ones(length(panel1),1),-(1-exp(-kappa1.* panel1))./kappa1];

    dtau2 = (r-alpha2+lambda2./kappa2 + sigma22.^2./(2.*kappa2.^2)-sigma21.*sigma22.*rho34./kappa2).*panel2+...
    sigma22.^2.*((1-exp(-2.*kappa2.*panel2))./(kappa2.^3) )./4+...
    (alpha2.*kappa2-lambda2+sigma21.*sigma22.*rho34-sigma22.^2./kappa2).*((1-exp(-kappa2.*panel2))./(kappa2.^2));
    
    Ztau2=[ones(length(panel2),1),-(1-exp(-kappa2.* panel2))./kappa2];

    dtau=[dtau1;dtau2];
    Ztau=[[Ztau1;zeros(length(panel2),2)],[zeros(length(panel1),2);Ztau2]];
end
function [logL,varargout]=kalmanFilter(P0,a0,c,T,varEta,d,Z,varEps,y)
%%
% THE TRANSITION EQUATION 
% a(t)=c+T*a(t-1)+eta(t)   eta~N(0,Q), Q=varEta
%
% THE MEASUREMENT EQUATION 
% y(t,tau)=d(tau)+Z(tau)*a(t)+epsilon(t)  epsilon~N(0,H), H=varEps


Ptt = P0;
att = a0; 
[nobs,N]=size(y);
save_att=zeros(size(y,1),length(a0));
save_dFtt_1=zeros(size(y,1),1);
save_vFv=zeros(size(y,1),1);
for t = 1:size(y,1)
    Ptt_1   = T*Ptt*T'+varEta;
    Ftt_1   = Z*Ptt_1*Z'+varEps;
    dFtt_1  = det(Ftt_1);
    invFtt_1 = inv(Ftt_1);

    att_1   = T*att + c;
    yt      = y(t,:)';
    ytt_1   = Z*att_1+d;
    vt      = yt-ytt_1;

    att = att_1 + Ptt_1*Z'*invFtt_1*(vt);
    Ptt = Ptt_1 - Ptt_1*Z'*invFtt_1*Z*Ptt_1;
    
    % ytt = Z*att+d;
    % vtt  = yt-ytt;

    % save_vtt(t,:) = vtt';
    % save_vt(t,:)    = (vt)';
    save_att(t,:)   = att';
    % save_Ptt_1(t,:) = [Ptt_1(1,1), Ptt_1(1,2), Ptt_1(2,1), Ptt_1(2,2)]; 
    % save_Ptt(t,:)   = [Ptt(1,1), Ptt(1,2), Ptt(2,1), Ptt(2,2)];

    save_dFtt_1(t,:)= dFtt_1;
    save_vFv(t,:)   = vt'*invFtt_1*vt;
    
end

    
logL = -(-(N*nobs/2)*log(2*pi)-0.5*sum(log(save_dFtt_1))-0.5*sum(save_vFv));

if nargout>1
    varargout{1}=save_att;
end
end
function [mu1,sigma11,sigma12,kappa1,alpha1,lambda1,mu2,sigma21,sigma22,kappa2,alpha2,lambda2,rho12,rho13,rho14,rho23,rho24,rho34,s,varargout]=paramUnpack(params)
    mu1=params(1);
    sigma11=params(2);
    sigma12=params(3);
    kappa1=params(4);
    alpha1=params(5);
    lambda1=params(6);
    mu2=params(7);
    sigma21=params(8);
    sigma22=params(9);
    kappa2=params(10);
    alpha2=params(11);
    lambda2=params(12);
    rho12=params(13);
    rho13=params(14);
    rho14=params(15);
    rho23=params(16);
    rho24=params(17);
    rho34=params(18);
    s=params(19:end);
    if nargout > 19
        paramDesc={'mu1','sigma11','sigma12','kappa1','alpha1','lambda1','mu2','sigma21','sigma22','kappa2','alpha2','lambda2','rho12','rho13','rho14','rho23','rho24','rho34'};
        for i=1:length(s)
            paramDesc{end+1}=sprintf('s%d',i);
        end
        varargout{1}=paramDesc;
    end
end