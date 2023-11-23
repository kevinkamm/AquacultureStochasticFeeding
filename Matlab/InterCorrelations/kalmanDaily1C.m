function [params,ss_att,negLogLikelihood,varargout]=kalmanDaily1C(dataName,dataStartDate,dataEndDate,dataPanel,r,varargin)

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


%% Extract data from spreadsheet
verb(@()disp('Load data'),verbose)
% dataName='Salmon';
data=load(['Data/Daily/',dataName]);
prices=data.output;
dates=data.udates;
ttm=data.uTTM;


dateInd = dates >=dataStartDate & dates <=dataEndDate;
ttmInd = ttm >= 0 & ttm <= dataPanel(end);

filledPrices=fillmissing(prices,'nearest',2);

yData=zeros(sum(dateInd),size(dataPanel,1));
panel=zeros(1,size(dataPanel,1));
for p=1:size(dataPanel,1)
    indP= ttm>dataPanel(p,1) & ttm<=dataPanel(p,2);
    yData(:,p)=mean(filledPrices(dateInd,indP),2);
    panel(p)=mean(ttm(indP));
end


yData=log(yData);
indContracts=1:size(yData,2); % start from 2 to exclude spot
% r=0.0303;

% dt = 1/365; % daily time step
dt = years(mean(diff(dates(dateInd)))); % =250 daily time step
% dt=1/300;

%%
disp('Calibrate')
verb(@()disp('Calibrate'),verbose)
% varEpsMode='single';
varEpsMode='diag';
% varEpsMode='tridiag';
ticEM=tic;
[params,negLogLikelihood]=EMschwartz(yData(:,indContracts),r,dt,panel(indContracts),varEpsMode);
ctimeEM=toc(ticEM);
verb(@()fprintf('Elapsed time for EM-Algo: %g s with log-likelihood %g\n',ctimeEM,-negLogLikelihood),verbose)

[P0,a0,c,T,varEta,dtau,Ztau,varEpsMode,~]=schwartzSSM(params,yData(:,indContracts),r,dt,panel(indContracts),varEpsMode);
[logL,ss_att]=kalmanFilter(P0,a0,c,T,varEta,dtau,Ztau,varEpsMode,yData(:,indContracts));

%%
if verbose > 0
    paramDesc=cell(length(params),1);
    paramDesc{1}='\mu';paramDesc{2}='\sigma_1';paramDesc{3}='\sigma_2';
    paramDesc{4}='\kappa';paramDesc{5}='\alpha';paramDesc{6}='\lambda';
    paramDesc{7}='\rho';
    for i=1:length(paramDesc)-7
        paramDesc{7+i}=sprintf('s%d',i);
    end
    paramTable=array2table(params,'VariableNames',paramDesc);
    disp(paramTable)
    
    %%
    if verbose >1
        fig=newFigure();
    
        plot(exp(yData(2:end,1)),'k','linewidth',1);
        plot(exp(yData(2:end,2:end)),'r--','linewidth',1);
        plot(exp(ss_att(:,1)),'b','linewidth',1);
    end
    %%
    % delta0=r is a good enough approx
    % disp('Soy parameters: mu, sigma1, sigma2, kappa, alpha, lambda, rho, delta0, Price0')
    disp(join(split(num2str([params(1:7),ss_att(end,2),exp(ss_att(end,1))]),whitespacePattern),', '))
    disp(join(split(num2str([params(1:7),ss_att(1,2),exp(ss_att(1,1))]),whitespacePattern),', '))
end

if nargout>3
    varargout{1}=dates(dateInd);
end

end
function verb(func,flag)
    if flag>0
            func();
    end
end

function [params,negLogLikelihood]=EMschwartz(yData,r,dt,panel,varEpsMode)
    % fminsearch does not work
    % fmincon changes with different bounds
    % fminunc with exp(params)->positive and rho=tanh(params)\in [0,1]
    %   similar to fmincon with reasonable bounds. During optimization many
    %   warnings with rcond
    % lsqnonlin does not make sense here
    % params=[mu,sigma11,sigma12,kappa1,alpha1,lambda1,rho]

    [nObservations,nContracts] = size(yData);
    
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

    % lb=zeros(1,7+n);
    lb=.01.*ones(1,7+n);lb(8:end)=0.001;
    lb(5)=-5;%alpha
    lb(7)=-1;
    ub=5.*ones(1,7+n);
    ub(7:end)=1;%ub(5)=1;
    x0=.1.*ones(1,7+n);
    x0(7)=.5;

    optionsF = optimoptions('fmincon','Display','off','UseParallel',false,'StepTolerance',1e-12,'OptimalityTolerance',1e-9,'MaxFunctionEvaluations',1e6);
    optionsU = optimoptions('fminunc','Display','off','UseParallel',false,'StepTolerance',1e-12,'OptimalityTolerance',1e-9,'MaxFunctionEvaluations',1e6);
    % optionsGA = optimoptions('ga','Display','iter','UseParallel',true);
    
    % x0 = ga(@(params) objective(params,yData,r,dt,panel),length(lb),[],[],[],[],lb,ub./10,@(params)mycon(params,yData,r,dt,panel),[],optionsGA);
    % x0 = ga(@(params) objective(params,yData,r,dt,panel),length(lb),[],[],[],[],lb,ub,[],[],optionsGA);
    [params,negLogLikelihood]= fmincon(...
        @(params) objective(params,yData,r,dt,panel,varEpsMode),...
        x0,[],[],[],[],lb,ub,[],optionsF);
    % [params,negLogLikelihood]= fminunc(...
    %     @(params) objective(params,yData,r,dt,panel,varEpsMode),...
    %     x0,optionsU);

end
function logL=objective(params,yData,r,dt,panel,varEpsMode)
    [P0,a0,c,T,varEta,dtau,Ztau,varEps,yData2]=schwartzSSM(params,yData,r,dt,panel,varEpsMode);
    logL=kalmanFilter(P0,a0,c,T,varEta,dtau,Ztau,varEps,yData2);
end
function [logL,ceq] = mycon(params,yData,r,dt,panel,varEpsMode)
    logL=-objective(params,yData,r,dt,panel,varEpsMode);
    ceq=[];
end
function [P0,a0,c,T,varEta,dtau,Ztau,varEps,yData2]=schwartzSSM(params,yData,r,dt,panel,varEpsMode)
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
    panel=panel(:); % time to maturities for different contracts in panel
    % params=exp(params);params(7)=tanh(log(params(7)));
    mu=params(1);
    sigma1=params(2);
    sigma2=params(3);
    kappa=params(4);
    alpha=params(5);
    lambda=params(6);
    rho=params(7);
    d=length(panel);
    
    switch varEpsMode
        case 'single'
            varEps= params(8).^2.*eye(d);
        case 'diag'
            varEps= diag(params(8:end).^2);
        case 'tridiag'
            s=reshape(1:d^2,d,d);
            ind=triu(s);ind=ind(ind>0);
            s=zeros(d,d);
            s(ind)=params(8:end);
            varEps=s'*s;
    end


    [nObservations,nContracts] = size(yData);
    
    varEta = [sigma1^2,              rho * sigma1 * sigma2; 
              rho * sigma1 * sigma2,               sigma2^2].*dt;

    c=[(mu-sigma1^2/2);kappa * alpha] .* dt;
    T=[1, -dt; 0, 1- kappa * dt];

    dtau = (r-alpha+lambda./kappa + sigma2.^2./(2.*kappa.^2)-sigma1.*sigma2.*rho./kappa).*panel+...
    sigma2.^2.*((1-exp(-2.*kappa.*panel))./(kappa.^3) )./4+...
    (alpha.*kappa-lambda+sigma1.*sigma2.*rho-sigma2.^2./kappa).*((1-exp(-kappa.*panel))./(kappa.^2));
    
    Ztau=[ones(nContracts,1),-(1-exp(-kappa.* panel))./kappa];

    P0=zeros(2,2);
    % a0=[mean(yData(:,1),1);r];
    a0=[yData(1,1);r];
    yData2=yData(2:end,:);
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

save_dFtt_1(save_dFtt_1<0)=0;%prevent complex number for unreasonable parameters
logL = -(-(N*nobs/2)*log(2*pi)-0.5*sum(log(save_dFtt_1))-0.5*sum(save_vFv));

if nargout>1
    varargout{1}=save_att;
end
end