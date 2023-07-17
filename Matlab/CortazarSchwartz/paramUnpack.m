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