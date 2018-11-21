function [y,Y,P,Y1] = ut(X,Wm,Wc,n,Rn,opt,Ts,u)
%Unscented Transformation
%Input:
%        f: nonlinear map
%        X: sigma points
%       Wm: weights for mean
%       Wc: weights for covraiance
%        n: numer of outputs of f
%        R: additive covariance
%Output:
%        y: transformed mean
%        Y: transformed sampling points
%        P: transformed covariance
%       Y1: transformed deviations

L = size(X,2);
y = zeros(n,1);
Y = zeros(n,L);

if opt == 1
    for k = 1:L
%         disp('entrando/saindo')
%         (X(:,k))'
        Y(:,k) = statef_est(X(:,k),Ts,u);   
%         (Y(:,k))'
%         pause
        y = y + Wm(k)*Y(:,k);       
    end
elseif opt == 2
    for k = 1:L                   
        Y(:,k) = measf(X(:,k));       
        y = y + Wm(k)*Y(:,k);       
    end
end

Y1 = Y - y(:,ones(1,L));


P = Y1*diag(Wc)*Y1'+Rn;
% if opt==1
%     y'
%     X
%     Y
%     Y1
%     pause
% end
P = (P+P')/2;
% 

%  opt
%  display('X')
%  size(X)
%  display('Wm')
%  size(Wm)
%  display('Wc')
%  size(Wc)
%  display('n')
%  n
%  display('Rn')
%  size(Rn)
%  display('y')
%  size(y)
%  display('Y1')
%  size(Y1)
%  display('diag Wc')
%  size(diag(Wc))
%  pause

 