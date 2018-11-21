
function [x3,P3,A3] = ekf(x3,P3,z3,Q3,R3,Ts,u)
%[x, P] = ekf(x,P,z,Qn,Rn,Ts,u);
% EKF   Extended Kalman Filter for nonlinear dynamic systems
% [x, P] = ekf(f,x,P,h,z,Q,R) returns state estimate, x and state covariance, P 
% for nonlinear dynamic system:
%           x_k+1 = f(x_k) + w_k
%           z_k   = h(x_k) + v_k
% where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
%       v ~ N(0,R) meaning v is gaussian noise with covariance R
% Inputs:   f: function handle for f(x)
%           x: "a priori" state estimate
%           P: "a priori" estimated state covariance
%           h: fanction handle for h(x)
%           z: current measurement
%           Q: process noise covariance 
%           R: measurement noise covariance
% Output:   x: "a posteriori" state estimate
%           P: "a posteriori" state covariance
%
% Example:
% 
% n=3;      %number of state
% q=0.1;    %std of process 
% r=0.1;    %std of measurement
% Q=q^2*eye(n); % covariance of process
% R=r^2;        % covariance of measurement  
% f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))];  % nonlinear state equations
% h=@(x)x(1);                               % measurement equation
% s=[0;0;1];                                % initial state
% x=s+q*randn(3,1); %initial state          % initial state with noise
% P = eye(n);                               % initial state covraiance
% N=20;                                     % total dynamic steps
% xV = zeros(n,N);          %estmate        % allocate memory
% sV = zeros(n,N);          %actual
% zV = zeros(1,N);
% for k=1:N
%   z = h(s) + r*randn;                     % measurments
%   sV(:,k)= s;                             % save actual state
%   zV(k)  = z;                             % save measurment
%   [x, P] = ekf(f,x,P,h,z,Q,R);            % ekf 
%   xV(:,k) = x;                            % save estimate
%   s = f(s) + q*randn(3,1);                % update process 
% end
% for k=1:3                                 % plot results
%   subplot(3,1,k)
%   plot(1:N, sV(k,:), '-', 1:N, xV(k,:), '--')
% end
%}
% By Yi Cao at Cranfield University, 02/01/2008
%

[x13,A3] = jac(x3,1,Ts,u);    %nonlinear update and linearization at current state
P3 = A3*P3*A3'+Q3;               %partial update
[z13,H3] = jac(x13,2,Ts,0);   %nonlinear measurement and linearization
P12 = P3*H3';   

%cross covariance
%K3 = P12*inv(H3*P12+R3);    %Kalman filter gain
%x3 = x13+K3*(z3-z13);       %state estimate
%P3 = P3 - K3*P12';

%state covariance matrix
R = chol(H3*P12+R3);                %Cholesky factorization
U = P12/R;                        %K=U/R'; Faster because of back substitution
x3 = x13+U*(R'\(z3-z13));           %Back substitution to get state update
P3 = P3-U*U';                     %Covariance update, U*U'=P12/R/R'*P12'=K*P12.
end

