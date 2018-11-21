function f = statef(x,Ts,u)
%f = [x(1)+Ts*10*(x(2)-x(1));x(2)+Ts*(28*x(1)-x(2)-x(1)*x(3));x(3)+Ts*(x(1)*x(2)-(8/3)*x(3))]+Ts*u;
% parametros máquina síncrona PSIM
Rs = 32.5;
%Rs = x(8);
%Rfd = 0.016;
Rfd = 346.8;
%Rfd = x(6);
%Rkd = 0.17;
%Rfd = x(9);
%Rkq = 0.17;
%Rkq = x(10);
%Lls = x(10);

%Lmq = 0.002;
%Lmq = 0.221;
%Lmd = 0.221;
Lmd = .889;
Lmq = 0.611;
Lls = 0.6*Lmd;
%Lmq = x(6);
%Lmd = 0.002;
%Lmd = Lmq;
Llfd = 30.3;
% Llkq = 0.00028;
% Llkd = 0.00091;
J = 0.001;
%J = x(9);
%B = 0.1;
B = 0.0;
Poles = 4;

% Com enrolamentos de amortecimento
% L_dq = [-(Lls+Lmq) 0           Lmq         0           0        0   0;
%      0          -(Lls+Lmd)  0           Lmd         Lmd         0   0;
%      -Lmq       0           (Llkq+Lmq)  0           0           0   0;
%      0          -Lmd        0          (Llkd+Lmd)   Lmd         0   0;
%      0          -Lmd        0           Lmd         (Llfd+Lmd)  0   0;
%      0          0           0           0           0           1.   0;
%      0          0           0           0           0           0   1.];

% Sem enrolamentos de amortecimento
L_dq = [-(Lls+Lmq) 0        0        0   0;
     0          -(Lls+Lmd)  Lmd      0   0;
     0          -Lmd        (Llfd+Lmd)  0   0;
     0          0           0           1.   0;
     0          0           0           0   1.];
 

% Com enrolamento de amortecimento
% R_dq = [ Rs                        0.5*Poles*x(6)*(Lls+Lmd)  0                   -0.5*Poles*x(6)*Lmd   -0.5*Poles*x(6)*Lmd  0        0;
%          -0.5*Poles*x(6)*(Lls+Lmq)    Rs                      0.5*Poles*x(6)*Lmq    0                   0                  0        0;
%          0                          0                       -Rkq                0                   0                  0        0;
%          0                          0                       0                   -Rkd                0                  0        0;
%          0                          0                       0                   0                   -Rfd               0        0;
%          0                          0                       0                   0                   0                  -B/J       0;
%          0                          0                       0                   0                   0                  0.5*Poles    0];

% Sem enrolamento de amortecimento
R_dq = [ Rs                        0.5*Poles*x(4)*(Lls+Lmd)  -0.5*Poles*x(4)*Lmd  0        0;
        -0.5*Poles*x(4)*(Lls+Lmq)    Rs                        0                  0        0;
         0                          0                       -Rfd               0        0;
         0                          0                         0                  -B/J       0;
         0                          0                         0                  0.5*Poles    0];

% Com enrolamento de amortecimento
% Te =(3/2)*(Poles/2)*( Lmd*(-x(2)+x(4)+x(5))*x(1) - Lmq*(-x(1)+x(3))*x(2) );  

% Sem enrolamento de amortecimento
% Te =(3/2)*(Poles/2)*(Lmd*(-x(2)+x(3))*x(1) - Lmq*(-x(1))*x(2));

% theta_e = x(end);
% while((theta_e > 2*pi) || (theta_e < 0))
%         theta_e = theta_e - sign(theta_e)*2*pi;
% end

f = x + Ts*(L_dq\(R_dq*x + u));


 