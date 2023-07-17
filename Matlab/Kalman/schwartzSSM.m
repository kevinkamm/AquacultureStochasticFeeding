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