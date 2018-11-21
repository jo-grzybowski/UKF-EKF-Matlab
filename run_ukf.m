function [t,sV,xV,zV,P_cov] = run_ukf(t_Ts,U,S)
% número de estados estimados
n_e = input('Número de parâmetros a estimar: ');
disp(' ')
disp('(1) simulação (2) dados')
opt = input('Tipo de entrada:  ');
format long

%Ts = 0.0001;
Ts = t_Ts(2) - t_Ts(1);
N_Ts = length(t_Ts);
disp('Taxa de amostragem: ')
disp(1/Ts)

rel = 1; % relação Ts/h

if opt == 1
    h = Ts/rel;
elseif opt == 2;
    h = Ts;
end

n = 5;      % number of states
nh = 5;     % number of output variables
q = 0.03;   % standard deviation of process 
%q = 1.692429537730919e-04;
r = 0.01;

%r = 1;
%r = [6 0 0 0;0 6 0 0;0 0 1 0;0 0 0 .01*2*pi];    %standard deviation of measurement
%r = [6 0 0 0;0 6 0 0;0 0 1 0;0 0 0 .01*2*pi*60];    %standard deviation of measurement

%r = 1.692429537730919e-04;
Qn = eye(n+n_e);

Qn = q^2*Qn; % covariance of process
Qn(1,1) = 0.2;
Qn(2,2) = 0.06;
%Q(3,3) = 0.01;
Rn = r^2*eye(nh,nh);        % covariance of measurement  

if opt == 1
    % Tensão de terminal
    Vd = 0;
    Vq = 100;
    % tensão de campo
    Vf = 10;
end

% frequencia
f = 60;

% Number of poles
Poles = 4;


% Sem enrolamentos
s = [0;0;0;0;0];                               % initial state
x = zeros(n+n_e,1); %initial state          % initial state with noise
%if n_e>0
%    x(6,1) = .5; 
%    x(7,1) = .5; 
%end
%x(8,1) = .15; 
P = eye(n+n_e);
if opt==2
    t = t_Ts;
else
    t = t_Ts(1):h:t_Ts(end);
end

N = length(t);                         % total dynamic 
xV = zeros(n+n_e,N_Ts);          %estmate        % allocate memory

if opt == 1
    sV = zeros(n,N);          %actual
elseif opt == 2
    sV = S;
end

zV = zeros(nh,N_Ts);

k1 = 0; % contador interno do laço for para amostragem
P_cov = zeros(n+n_e,n+n_e,N);

%initializes the input vector
u = zeros(n+n_e,1);
%rQn = [.1 -.1 .1 -.1 .1 -.1 .1]';
%rRn = [.1 -.1 .1 -.1]';
for k = 1:N
    % update input vectorqn
     if opt == 1
        u(1) = Vq;
        u(2) = Vd;
        u(3) = Vf;
        s = statef(s,h,u(1:n));
        s = s + q*randn(n,1);      % update process 
     end

    if mod(h*(k-1),Ts) == 0
        k1 = k1 + 1;
        if opt == 2
            % com enrolamentos
%             u(1) = U(1,k1);
%             u(2) = U(2,k1);
%             u(5) = U(3,k1);
%             s(1) = S(1,k);
%             s(2) = S(2,k);
%             s(5) = S(5,k);
%             s(7) = S(7,k);

            % sem enrolamentos
            u(1) = U(1,k1);
            u(2) = U(2,k1);
            u(3) = U(3,k1);
            s(1) = S(1,k);
            s(2) = S(2,k);
            s(3) = S(3,k);
            s(4) = S(4,k)/2; % sobre dois porque é a w do rotor
            s(5) = S(5,k);
        end
        z = measf(s);
        if opt == 1
            z = z + r*randn(nh,1);              % measurments                            % save actual state
        end

        zV(:,k1)  = z;                           % save measurment
        [x, P] = ukf(x,P,z,Qn,Rn,Ts,u);           % ukf 
        P_cov(:,:,k1) = P;
        xV(:,k1) = x;  % save estimate
    end
    sV(:,k) = s;
end

% Plot State variables
save statevars_ukf
% com enrolamentos de amortecimento
%plotar(t,t_Ts,xV,sV,zV,P_cov)

% sem enrolamentos de amortecimento

plotar5estados(t,t_Ts,xV,sV,zV,P_cov)
  