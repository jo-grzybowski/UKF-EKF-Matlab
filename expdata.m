function [t,U,S,If] = expdata(sample_freq)
% SC1: data from scope #1
% SC2: data from scope #2
filename1 = 't0059.txt';
filename2 = 't0060.txt';
SC1 = importdata(filename1);
SC2 = importdata(filename2);
dados = [SC1.data SC2.data];

% características construtivas
polos = 4;

% colunas da matriz 'dados' original:
% 1 t
% 2 VA
% 3 VA_PICO
% 4 VB
% 5 VB_PICO
% 6 VC
% 7 VC_PICO
% 8 THETA
% 9 THETA_PICO
% 10 t
% 11 IA
% 12 IA_PICO
% 13 IB
% 14 IB_PICO
% 15 IC
% 16 IC_PICO

% eliminar colunas
elim = [20 18 16 14 12 11 10 9 7 5 3];
for k = 1:length(elim)
    dados(:,elim(k)) = [];
end

% colunas da matriz 'dados' nova:
% 1 t
% 2 VA
% 3 VB
% 4 VC
% 5 THETA
% 6 IA
% 7 IB
% 8 IC
% 9 IF

% elimina offset das colunas de tensão AC e corrente AC
col = [2 3 4 6 7 8];
for k = 1:length(col)
    dados(:,col(k)) = dados(:,col(k));
end

% ajuste da frequencia de amostragem
freq = 1/(dados(2,1)-dados(1,1));
sampler = round(freq/sample_freq);
N1 = size(dados,1);
N2 = size(dados,2);
dados_sampled = zeros(round(N1/sampler)+1,N2);
k1 = 0; % counter
for k = 1:sampler:N1
    k1 = k1+1;
    dados_sampled(k1,:) = dados(k,:);
end
dados_sampled(end,:) = dados(end,:);

% Lfd = T0011.txt, 13V
% Ls = T0012.txt, 13V
% dados amostrados na frequencia solicitada
% atribuição das variáveis conforme montagem
fator_I = (1/0.100)/10; % Alicate amperímetro configurada para 100mV por Ampère e ponteira osciloscópio 10x.
fator_V = (1/0.01)/10; % Voltímetro setado para 100:1 e osciloscópio
dc_offset = 0.16; % ajuste no OFFSET DC da ponteira de corrente
t = dados_sampled(:,1);
Va = fator_V*dados_sampled(:,4);
Vb = fator_V*dados_sampled(:,2);
Vc = fator_V*dados_sampled(:,3);

Ia = fator_I*dados_sampled(:,8);
Ib = fator_I*dados_sampled(:,6);
Ic = fator_I*dados_sampled(:,7);
If = fator_I*dados_sampled(:,9);
lead = 20; % moving average filter lead
lag = 20;  % moving average filter lag
If = movavg(If,lead,lag);

If = If - dc_offset; 

theta_m = dados_sampled(:,5);

% trabalhar com os dados de ângulo
lim1 = max(theta_m)*0.32;
theta_mf = theta_m > lim1;
theta_s1 = zeros(length(theta_mf),1);
sweep1 = 10;

for k = sweep1:length(theta_mf)-sweep1
    theta_s1(k) = sum(theta_mf(k-sweep1+1:k+sweep1));
end
lim2 = max(theta_s1)*0.7;
theta_s1 = theta_s1 > lim2;
theta_s2 = zeros(length(theta_mf),1);
sweep2 = 5;

for k = sweep2:length(theta_mf)-sweep2
    theta_s2(k) = sum(theta_s1(k-sweep2+1:k+sweep2));
end
lim3 = max(theta_s2)*0.7;
theta_s2 = theta_s2 > lim3;

theta_s3 = zeros(length(theta_mf),1);
sweep3 = 5;

for k = sweep3:length(theta_mf)-sweep3
    theta_s3(k) = sum(theta_s1(k-sweep3+1:k+sweep3));
end
lim4 = max(theta_s3)*0.9;
theta_s3 = theta_s3 > lim4;

P = Pfind(find(theta_s3));
theta_pulse = zeros(length(theta_s3),1);
for k = 1:length(P)
    theta_pulse(P(k)) = 1;
end

% theta_e --> ângulo elétrico
%orient = 1; % negativa
orient = -1; %positiva
theta_e = zeros(length(theta_mf),1);
theta_m = zeros(length(theta_mf),1);
for k = 1:length(P)-1
    theta_e(P(k):P(k+1)) = linspace(0,orient*(-(polos/2)*(2*pi)),P(k+1)-P(k)+1);
    theta_m(P(k):P(k+1)) = linspace(0,orient*(-(2*pi)),P(k+1)-P(k)+1);
end

% agora antes do primeiro e após o último sinal de pulso

step = abs(theta_e(P(1))-theta_e(P(1)+1));
step_m = abs(theta_m(P(1))-theta_m(P(1)+1));
for k = P(1):-1:1
    theta_e(k) = orient*(-(polos/2)*(2*pi)+step*(P(1) - k));
    theta_m(k) = orient*(-(2*pi)+step_m*(P(1) - k));
end

step = abs(theta_e(P(end))-theta_e(P(end)-1));
step_m = abs(theta_m(P(end))-theta_m(P(end)-1));
for k = P(end):length(theta_e)
    theta_e(k) = orient*((-(polos/2)*(2*pi)-step*(k - P(end))));
    theta_m(k) = orient*((-(2*pi)-step_m*(k - P(end))));
end
theta_e = mod(theta_e,2*pi);
theta_m = mod(theta_m,2*pi);

theta_eint = zeros(1,length(theta_e));
theta_mint = zeros(1,length(theta_m));

% theta continuo para filtro
theta_eint(1) = theta_e(1);
theta_mint(1) = theta_m(1);

if (theta_pulse(2)~=1 & theta_pulse(1)~=1)
    step = theta_e(2) - theta_e(1);
    step1 = theta_m(2) - theta_m(1);
end

for k = 2:length(theta_e);
    if (theta_pulse(2)~=1 & theta_pulse(1)~=1)
        step = theta_e(2) - theta_e(1);
        step1 = theta_m(2) - theta_m(1);
    end
    theta_eint(k) = theta_eint(k-1) + step;
    theta_mint(k) = theta_mint(k-1) + step1;
end

% velocidade angular para os gráficos
w = zeros(1,length(theta_e));
dt = t(2)-t(1);
for k = 2:length(theta_e)-1
    w(k)= (theta_eint(k+1) - theta_eint(k-1))/(2*dt);
end
w(1) = w(2);
w(end) = w(end-1);

% Transformacao ABC-DQ0
Vdqn = dqn([Va';Vb';Vc'],theta_e);
Vd = Vdqn(1,:);
Vq = Vdqn(2,:);
V0 = Vdqn(3,:);

Idqn = dqn([Ia';Ib';Ic'],theta_e);
Id = Idqn(1,:);
Iq = Idqn(2,:);
I0 = Idqn(3,:);


% Figures
figure(1)
set(gcf,'Position',[0 0 850 400])
subplot(2,1,1)
plot(t,Va,t,Vb,t,Vc)
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('Va','Vb','Vc')
grid on

subplot(2,1,2)
plot(t,Ia,t,Ib,t,Ic)
xlabel('Time (s)')
ylabel('Current (A)')
legend('Ia','Ib','Ic')
grid on

figure(2)
set(gcf,'Position',[0 0 850 400])
subplot(2,1,1)
[AX,H1,H2] = plotyy(t,theta_m,t,theta_pulse)

xlabel('Time (s)')
set(get(AX(1),'Ylabel'),'String','Rotor q-axis position')
set(get(AX(2),'Ylabel'),'String','Refined position pulse')
set(H1,'LineWidth',2)
set(H2,'LineWidth',2)
set(H2,'LineStyle','--')
grid on

subplot(2,1,2)
[AX,H1,H2] = plotyy(t,theta_pulse,t,theta_e)
xlabel('Time (s)')
set(get(AX(1),'Ylabel'),'String','Rotor q-axis position')
set(get(AX(2),'Ylabel'),'String','\theta_e')
set(H1,'LineWidth',2)
set(H2,'LineWidth',2)
set(H2,'LineStyle','--')
grid on

figure(3)
set(gcf,'Position',[0 0 850 400])

subplot(2,1,1)
plot(t,Va,t,Vb,t,Vc)
legend('Va','Vb','Vc')
xlabel('Time (s)')
ylabel('Voltage (V)')

subplot(2,1,2)
plot(t,Vd,t,Vq,t,V0)
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('Vd','Vq','V0')

figure(4)
set(gcf,'Position',[0 0 850 400])
subplot(2,1,1)
plot(t,Ia,t,Ib,t,Ic)
legend('Ia','Ib','Ic')
xlabel('Time (s)')
ylabel('Current (A)')

subplot(2,1,2)
plot(t,Id,t,Iq,t,I0)
legend('Id','Iq','I0')
xlabel('Time (s)')
ylabel('Current (A)')

figure(5)
set(gcf,'Position',[0 0 600 400])
plot(t,If)
xlabel('Time (s)')
ylabel('Field current (A)')
ylim([0 0.5])
grid on

figure(6)
plot(t,theta_e,t,theta_m)
xlabel('Time (s)')
ylabel('Angular position of the rotor (A)')

% variáveis de saída para entrada na UKF, EKF
Vf = 64;
Rfd = 346.8;

S = zeros(5,length(t));
S(1,:) = Iq;
S(2,:) = Id;
%S(3,:) = Vf/Rfd*ones(1,length(t));
S(3,:) = If;
S(4,:) = w;
S(5,:) = theta_mint;

U = zeros(3,length(t));
U(1,:) = Vq;
U(2,:) = Vd;
U(3,:) = Vf*ones(1,length(t));
end

function P = Pfind(P0)
n = length(P0);
P1 = zeros(n,1);
k = 1;
P1(1) = 1;
for j = 1:n-1
    if P0(j+1)-P0(j)==1
        P1(j+1) = k;
    else
        k = k+1;
        P1(j+1) = k;
    end
end
k1 = 1;
for j = 1:k
    vec = find(P1==j);
    P(k1) = floor(mean(P0(vec)));
    k1 = k1 + 1;
end

end

function Vdqn = dqn(Vabc,theta)
n = size(Vabc,2);
Vdqn = zeros(3,n);
for k = 1:n
    Vdqn(:,k) = (2/3)*[sin(theta(k)) sin(theta(k)-2/3*pi) sin(theta(k)+2/3*pi);cos(theta(k)) cos(theta(k)-2/3*pi) cos(theta(k)+2/3*pi);sqrt(2)/2 sqrt(2)/2 sqrt(2)/2]*Vabc(:,k);
end
end





