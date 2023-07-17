function [W1,W2,varargout]=brownianMotions(T,N,M,rho)
%%BROWNIANMOTIONS simulates M trajectories of two correlated Brownian
% motions W_t^1 and W_t^2 with correlation \rho in the interval [0,T] with
% N points in a homogeneous grid.
%
%   Input:
%
%   Output:
%   
%   Usage:
%       [W1,W2]=brownianMotions(T,N,M,rho)
%       [W1,W2,t]=brownianMotions(T,N,M,rho): additionally outputs time
%                                             grid t
%   See also:
%       randn, cumsum

    dt=T./(N-1);
    W=zeros(N,M,2);
    W(2:end,:,:)=sqrt(dt).*randn(N-1,M,2);
    W=cumsum(W,1);

    W1=W(:,:,1);
    W2=rho.*W1+sqrt(1-rho.^2).*W(:,:,2);
    
    if nargout>2
        varargout{1}=linspace(0,T,N)';
    end
end