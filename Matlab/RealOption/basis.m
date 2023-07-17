function a=basis(x,y)
    %% Polynomial
    %   x (m x nd double): nd number of points, m number of variables
    %   a (q x nd double): vector of all x_i*x_j for multivariate
    %                      polynomial with degree 2 evaluated at x
    %                      order: 1,x_i,x_i^2,x_i*x_j
    %                      q = 1+2*m+(m over 2) = const+lin+quad+mixed
    %                      if y is given, then the least square approx is 
    %                      returned

    sz=size(x);
    m=sz(1);
    nd=sz(2:end);
    n=prod(nd);
    C = nchoosek(1:m,2);% each row has indices for combination x_i*x_j
    a=ones([1+2*m+size(C,1),n]);
    x=reshape(x,m,n);
    a(2:m+1,:)=x;
    a(m+2:2*m+1,:)=x.^2;
    a(2*m+2:end,:)=x(C(:,1),:).*x(C(:,2),:);
    a=a';
    if nargin>1
        coeff=a'*a\(a'*y(:));%least square approx: min norm(basis(x)*(basis(x)\y)-y)
        a=a*coeff;
    end
    a=reshape(a,size(y));
end
% function a=basis(x,y)
%     %% Laguerre
%     %   x (m x nd double): nd number of points, m number of variables
%     %   a (q x nd double): vector of all x_i*x_j for multivariate
%     %                      polynomial with degree 2 evaluated at x
%     %                      order: 1,x_i,x_i^2,x_i*x_j
%     %                      q = 1+2*m+(m over 2) = const+lin+quad+mixed
%     %                      if y is given, then the least square approx is 
%     %  syms x;laguerreL(n,x)
%     f1=@(x) 1 - x;
%     f2=@(x) x.^2./2 - 2.*x + 1;
% 
%     sz=size(x);
%     m=sz(1);
%     nd=sz(2:end);
%     n=prod(nd);
%     C = nchoosek(1:m,2);% each row has indices for combination x_i*x_j
%     a=ones([1+2*m+size(C,1),n]);
%     x=reshape(x,m,n);
%     a(2:m+1,:)=f1(x);
%     a(m+2:2*m+1,:)=f2(x);
%     a(2*m+2:end,:)=f1(x(C(:,1),:)).*f1(x(C(:,2),:));
%     a=a';
%     if nargin>1
%         coeff=a'*a\(a'*y(:));%least square approx: min norm(basis(x)*(basis(x)\y)-y)
%         a=a*coeff;
%     end
%     a=reshape(a,size(y));
% end