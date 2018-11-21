function f = statef_est(x,Ts,u)
%f = [x(1)+Ts*10*(x(2)-x(1));x(2)+Ts*(28*x(1)-x(2)-x(1)*x(3));x(3)+Ts*(x(1)*x(2)-(8/3)*x(3))]+Ts*u;
% parametros máquina síncrona PSIM
%Rs = 0.1;
Rs = 32.5;
%Rs = x(6);
%Rfd = 0.016;
Rfd = 346.8;
%Rfd = x(7);
Rkd = 0.0;
%Rfd = x(9);
Rkq = 0.0;
%Rkq = x(10);
%Lls = x(10);

%Lmq = .22100;
%Lmd = .22100;
Lmd = x(6);
Lmq = x(7);
%Lmd = 0.73;
%Lmq = 0.58;
%Lmq = 0.7*Lmd;
%Lmd = 0.893;
%Lmq = 0.611;
%Lmq = x(6);
%Lls = x(6);
Lls = 0.15*Lmd;
%Lls = 0.15*Lmd;
%Lmd = 0.018;
%Lmq = x(6);
%Lmd = Lmq;
%Lmd = 0.002;

Llfd = 2.23;
%Llfd = x(6);
%Llfd = x(7);
Llkq = 0.0;
Llkd = 0.0;
J = 0.001;
%J = x(9);
%B = 0.1;
B = 0.001;
Poles = 4;

%J = x(9);
n_e = length(x) - 5;
                                                                        
L_dq = [-(Lls+Lmq) 0        0        0       0;
     0          -(Lls+Lmd)  Lmd         0       0;
     0          -Lmd        (Llfd+Lmd)  0       0;
     0          0           0           1       0;
     0          0           0           0       1];

% L_dq_inv = [- 1/(Lls + Lmq) - Lmq^2/((Lls + Lmq)^2*(Llkq + Lmq - Lmq^2/(Lls + Lmq))),0, Lmq/((Lls + Lmq)*(Llkq + Lmq - Lmq^2/(Lls + Lmq))),0,0, 0, 0;
%             0,- (Lmd/(Lls + Lmd) - (Lmd*(Lmd - Lmd^2/(Lls + Lmd)))/((Lls + Lmd)*(Llkd + Lmd - Lmd^2/(Lls + Lmd))))^2/(Llfd + Lmd - Lmd^2/(Lls + Lmd) - (Lmd - Lmd^2/(Lls + Lmd))^2/(Llkd + Lmd - Lmd^2/(Lls + Lmd))) - 1/(Lls + Lmd) - Lmd^2/((Lls + Lmd)^2*(Llkd + Lmd - Lmd^2/(Lls + Lmd))),                                                  0, Lmd/((Lls + Lmd)*(Llkd + Lmd - Lmd^2/(Lls + Lmd))) - ((Lmd - Lmd^2/(Lls + Lmd))*(Lmd/(Lls + Lmd) - (Lmd*(Lmd - Lmd^2/(Lls + Lmd)))/((Lls + Lmd)*(Llkd + Lmd - Lmd^2/(Lls + Lmd)))))/((Llkd + Lmd - Lmd^2/(Lls + Lmd))*(Llfd + Lmd - Lmd^2/(Lls + Lmd) - (Lmd - Lmd^2/(Lls + Lmd))^2/(Llkd + Lmd - Lmd^2/(Lls + Lmd)))), (Lmd/(Lls + Lmd) - (Lmd*(Lmd - Lmd^2/(Lls + Lmd)))/((Lls + Lmd)*(Llkd + Lmd - Lmd^2/(Lls + Lmd))))/(Llfd + Lmd - Lmd^2/(Lls + Lmd) - (Lmd - Lmd^2/(Lls + Lmd))^2/(Llkd + Lmd - Lmd^2/(Lls + Lmd))), 0, 0;
%             -Lmq/((Lls + Lmq)*(Llkq + Lmq - Lmq^2/(Lls + Lmq))),0,1/(Llkq + Lmq - Lmq^2/(Lls + Lmq)),0,0, 0, 0;
%             0, ((Lmd - Lmd^2/(Lls + Lmd))*(Lmd/(Lls + Lmd) - (Lmd*(Lmd - Lmd^2/(Lls + Lmd)))/((Lls + Lmd)*(Llkd + Lmd - Lmd^2/(Lls + Lmd)))))/((Llkd + Lmd - Lmd^2/(Lls + Lmd))*(Llfd + Lmd - Lmd^2/(Lls + Lmd) - (Lmd - Lmd^2/(Lls + Lmd))^2/(Llkd + Lmd - Lmd^2/(Lls + Lmd)))) - Lmd/((Lls + Lmd)*(Llkd + Lmd - Lmd^2/(Lls + Lmd))),                                                  0,                                                                                                                  1/(Llkd + Lmd - Lmd^2/(Lls + Lmd)) + (Lmd - Lmd^2/(Lls + Lmd))^2/((Llkd + Lmd - Lmd^2/(Lls + Lmd))^2*(Llfd + Lmd - Lmd^2/(Lls + Lmd) - (Lmd - Lmd^2/(Lls + Lmd))^2/(Llkd + Lmd - Lmd^2/(Lls + Lmd)))),                                      -(Lmd - Lmd^2/(Lls + Lmd))/((Llkd + Lmd - Lmd^2/(Lls + Lmd))*(Llfd + Lmd - Lmd^2/(Lls + Lmd) - (Lmd - Lmd^2/(Lls + Lmd))^2/(Llkd + Lmd - Lmd^2/(Lls + Lmd)))), 0, 0;
%             0,-(Lmd/(Lls + Lmd) - (Lmd*(Lmd - Lmd^2/(Lls + Lmd)))/((Lls + Lmd)*(Llkd + Lmd - Lmd^2/(Lls + Lmd))))/(Llfd + Lmd - Lmd^2/(Lls + Lmd) - (Lmd - Lmd^2/(Lls + Lmd))^2/(Llkd + Lmd - Lmd^2/(Lls + Lmd))),                                                  0,                                                                                                                                                          -(Lmd - Lmd^2/(Lls + Lmd))/((Llkd + Lmd - Lmd^2/(Lls + Lmd))*(Llfd + Lmd - Lmd^2/(Lls + Lmd) - (Lmd - Lmd^2/(Lls + Lmd))^2/(Llkd + Lmd - Lmd^2/(Lls + Lmd)))),                                                                                                  1/(Llfd + Lmd - Lmd^2/(Lls + Lmd) - (Lmd - Lmd^2/(Lls + Lmd))^2/(Llkd + Lmd - Lmd^2/(Lls + Lmd))), 0, 0;
%             0,0,0,0,0,1,0;
%             0,0,0,0,0,0,1];
      

L_dq = [L_dq zeros(5,n_e);zeros(n_e,5) eye(n_e,n_e)];

R_dq = [ Rs                        0.5*Poles*x(4)*(Lls+Lmd) -0.5*Poles*x(4)*Lmd     0           0;
         -0.5*Poles*x(4)*(Lls+Lmq)  Rs                        0                     0           0;
         0                          0                       -Rfd                    0           0;
         0                          0                       0                      -B/J        0;
         0                          0                       0                       0.5*Poles   0];
     
%        disp('R_dq')
%         R_dq
%         pause

R_dq = [R_dq zeros(5,n_e);zeros(n_e,5) eye(n_e)];

% Com enrolamentos de amortecimento
%Te =(3/2)*(Poles/2)*(Lmd*(-x(2)+x(4)+x(5))*x(1) - Lmq*(-x(1)+x(3))*x(2));  

% Sem enrolamentos de amortecimento
% Te =(3/2)*(Poles/2)*(Lmd*(-x(2)+x(3))*x(1) - Lmq*(-x(1))*x(2)); 

% theta_e = x(7);
% while((theta_e > 2*pi) || (theta_e < 0))
%         theta_e = theta_e - sign(theta_e)*2*pi;
% end

f = x + Ts*(L_dq\(R_dq*x + u));



