function [mu1,sigma11,sigma12,kappa1,alpha1,lambda1,rho,s,varargout]=paramUnpack(params)
    mu1=params(1);
    sigma11=params(2);
    sigma12=params(3);
    kappa1=params(4);
    alpha1=params(5);
    lambda1=params(6);
    rho=params(7);
    s=params(8:end);
    if nargout > 19
        paramDesc={'mu','sigma1','sigma2','kappa','alpha','lambda','rho'};
        for i=1:length(s)
            paramDesc{end+1}=sprintf('s%d',i);
        end
        varargout{1}=paramDesc;
    end
end