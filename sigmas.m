function X=sigmas(x,P,c)
%Sigma points around reference point
%Inputs:
%       x: reference point
%       P: covariance
%       c: coefficient
%Output:
%       X: Sigma points
P = (P+P')/2;


% flag1 = 0;
 while min(real(eig(P)))<0
     P = P + 1*eye(size(P,1),size(P,2));
%      disp('Eigenvalue...')
%      flag1 = 1;
 end
%  
%  if flag1 == 1
%      disp('Ok')
%      flag1 = 0;
%  end


A = c*chol(P)';
Y = x(:,ones(1,numel(x)));
X = [x Y+A Y-A];
% display('input sigma')
% P
% 
% c*chol(P)'
% pause


