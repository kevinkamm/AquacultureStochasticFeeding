function [A,B,C,D,Mean0,Cov0,StateType]=schwartzSSM(params,dt,panel)
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
%   State equation: [x(t+1),delta(t+1)] = c + T [x(t),delta(t)] + eta
%   Observation equation: y(t+1,tau) = d(tau) + Z(tau) y(t,tau) + epsilon
%   Var(eta) = [sigma_1^2,              rho sigma_1 sigma_2; 
%               rho sigma_1 sigma_2,               sigma_2^2]*dt
%   c=[(mu-sigma_1^2/2);kappa alpha] * dt
%   T = [1, -dt; 0, 1- kappa dt]
%   d(tau) = A(tau)
%   Z(tau) = [1,-(1-e^(-kappa tau))/kappa]
%
%
% Matlab SSM:
%   State equation:       x(t) = A * x(t-1) + B * u(t)
%   Observation equation: y(t) = C * x(t)   + D * e(t)
%
% To use Matlab built-in functions, add a state xi for the constants
%   x(t) -> [x(t),delta(t),xi(t)] with StateType=[0,0,1]
%   A -> [T, c; 0 1]
%   B -> [chol(Var(eta));0]
%   C -> [Z(tau), d(tau)]
%   D -> D

    StateType=[0,0,1];
    Mean0=[];
    Cov0=[];

    mu=params(1);
    sigma1=params(2);
    sigma2=params(3);
    kappa=params(4);
    alpha=params(5);
    lambda=params(6);
    rho=params(7);
    r=params(8);

    szPanel=length(panel);
    
    % Define A for Matlab
    A=zeros(3,3);
    A(1,1)=1;A(1,2)=-dt;A(1,3)=(mu-sigma1^2/2)*dt;
    A(2,2)=1-kappa*dt;A(2,3)=kappa*alpha*dt;
    A(3,3)=1;
    
    % Define B for Matlab
    % cholesky decomposition of var(eta)=stdEta*stdEta'
    stdEta=[sigma1,0;rho*sigma2, sigma2*sqrt(1-rho^2)];
    B=zeros(3,2);
    B(1:2,1:2)=stdEta;
    
    % Define C for Matlab
    dtau = (r-alpha+lambda./kappa + sigma2.^2./(2.*kappa.^2)-sigma1.*sigma2.*rho./kappa).*panel+...
    sigma2.^2.*((1-exp(-2.*kappa.*panel))./(kappa.^3) )./4+...
    (alpha.*kappa-lambda+sigma1.*sigma2.*rho-sigma2.^2./kappa).*((1-exp(-kappa.*panel))./(kappa.^2));
    
    Ztau=[ones(szPanel,1),-(1-exp(-kappa.* panel))./kappa];
    C=[Ztau,dtau];
    

    % Define D for Matlab
    D=[];%[NaN;NaN];
end