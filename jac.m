function [z,A] = jac(x4,opt,Ts,u)
% JACCSD Jacobian through complex step differentiation
% [z J] = jaccsd(f,x)
% z = f(x)
% J = f'(x)
%
if opt == 1
    z = statef_est(x4,Ts,u);
elseif opt ==2
    z = measf(x4);
end

n = length(x4);
m = length(z);
A = zeros(m,n);

incr = n*eps;
for k = 1:n
    x1 = x4;
    x1(k)  = x1(k)+incr*i;
    if opt == 1
        x = statef_est(x1,Ts,u);
    elseif opt == 2
        x = measf(x1);
    end
    A(:,k) = imag(x)/incr;
end
end