function [x,P] = ukf(x,P,z,Q,R,Ts,u)
% UKF   Unscented Kalman Filter for nonlinear dynamic systems
% [x, P] = ukf(f,x,P,h,z,Q,R) returns state estimate, x and state covariance, P 
% for nonlinear dynamic system (for simplicity, noises are assumed as additive):
%           x_k+1 = f(x_k) + w_k
%           z_k   = h(x_k) + v_k
% where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
%       v ~ N(0,R) meaning v is gaussian noise with covariance R
% Inputs:   f: function handle for f(x)
%           x: "a priori" state estimate
%           P: "a priori" estimated state covariance
%           h: function handle for h(x)
%           z: current measurement
%           Q: process noise covariance 
%           R: measurement noise covariance
% Output:   x: "a posteriori" state estimate
%           P: "a posteriori" state covariance
%
% Example:
%{
n=3;      %number of state
q=0.1;    %standard deviation of process 
r=0.1;    %standard deviattion of measurement
Q=q^2*eye(n); % covariance of process
R=r^2;        % covariance of measurement  
f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))];  % nonlinear state equations
h=@(x)x(1);                               % measurement equation
s=[0;0;1];                                % initial state
x=s+q*randn(3,1); %initial state          % initial state with noise
P = eye(n);                               % initial state covariance
N=20;                                     % total dynamic steps
xV = zeros(n,N);          %estmate        % allocate memory
sV = zeros(n,N);          %actual
zV = zeros(1,N);
for k=1:N
  z = h(s) + r*randn;                     % measurments
  sV(:,k)= s;                             % save actual state
  zV(k)  = z;                             % save measurment
  [x, P] = ukf(f,x,P,h,z,Q,R);            % ekf 
  xV(:,k) = x;                            % save estimate
  s = f(s) + q*randn(3,1);                % update process 
end
for k=1:3                                 % plot results
  subplot(3,1,k)
  plot(1:N, sV(k,:), '-', 1:N, xV(k,:), '--')
end
%}
%
% By Yi Cao at Cranfield University, 04/01/2008
%

L = length(x);                                %number of states
m = length(z);                                %number of measurements
alpha = 1e-3;                                 %default, tunable
ki = 0;                                       %default, tunable
beta = 2;                                     %default, tunable
lambda = alpha^2*(L+ki)-L;  
%scaling factor
c = L+lambda;  %scaling factor
Wm = [lambda/c 0.5/c+zeros(1,2*L)];           %weights for means
Wc = Wm;
Wc(1) = Wc(1)+(1-alpha^2+beta);               %weights for covariance
c = sqrt(c);
X = sigmas(x,P,c);                            %sigma points around x
% X
% pause
% clc
[x1,X1,Px,X2] = ut(X,Wm,Wc,L,Q,1,Ts,u);       %unscented transformation of process
% x1'
% pause
% clc
% X1
% pause
% clc
% X2
% pause
% clc
% Px
% pause
% clc
% pause
[z1,Z1,Py,Z2] = ut(X1,Wm,Wc,m,R,2,Ts,0);      %unscented transformation of measurments

%  z1'
%  pause
%  clc
%  Z1
%  pause
%  clc
%   Z2
%  pause
%  clc
%  Py
%  pause
%  clc
 
Pxy = X2*diag(Wc)*Z2'; %transformed cross-covariance
% Pxy
% pause
% clc
% inv(Py)
%  pause
%  clc
%R = chol(Py);
%K = (Pxy/R)/R';                               % Filter gain.
%Py
%inv(Py)
%pause
K = Pxy*inv(Py);
% K
% pause
% clc
% Pxy
% K

x = x1 + K*(z-z1);
% z'
% (z-z1)'
% K*(z-z1)
% pause
% x'
% pause
% clc

P = Px - K*(Py*K'); %covariance update
%  z'
%  x'
%  P
%  pause
%  clc




