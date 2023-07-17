function [P,delta]=schwartzTwoFactor(delta0,P0,r,sigma1,sigma2,kappa,alpha,lambda,W1,W2,t)
%%Simulates the Schwartz' two-factor model for given correlated Brownian
% motions W_t^1 and W_t^2 and parameters:
%   dP_t = (r-\delta_t) P_t dt + \sigma_1 P_t dW^1_t
%   d\delta_t = (\kappa (\alpha-\delta_t)-\lambda) dt + \sigma_2 dW^2_t
%   d<W^1,W^2>_t = \rho dt
%
%   Input:
%
%   Output:
%
%   Usage:
%       [P,delta]=schwartzTwoFactor(mu,sigma1,sigma2,kappa,alpha,W1,W2,t)
%
%   See also:
%       brownianMotions
%
    N=length(t);
    dt=t(end)./(N-1);
    dW2=diff(W2,1);

    % Ornstein-Uhlenbeck
    expkt=exp(-kappa.*t);
    delta=delta0.*expkt+...
            (alpha-lambda./kappa).*(1-expkt)+...
            expkt.*sigma2.*stochInt(exp(kappa.*t),dW2);

    % Geometrical Brownian motion
    Ideltadt=lebInt(delta,dt,N);
    P=P0.*exp((r-sigma1.^2./2).*t-Ideltadt+sigma1.*W1);

end
function Idt=lebInt(f,dt,N)
    [~,m]=size(f);
    Idt=zeros(N,m);
    Idt(2:end,:)=cumsum(f(1:end-1,:).*dt,1);
end
function IdW=stochInt(f,dW)
    [n,m]=size(dW);
    IdW=zeros(n+1,m);
    IdW(2:end,:)=cumsum(f(1:end-1,:).*dW,1);
end