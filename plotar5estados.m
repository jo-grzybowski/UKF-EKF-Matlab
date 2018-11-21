function plotar5estados(t,t_Ts,xV,sV,zV,P_cov)
 nome = input('Nome da simulação: ','s');
 
% Plot State variables
  figure(1)
  subplot(3,1,1)
  plot(t, sV(1,:), '-',t_Ts, xV(1,:), '--');
  set(gca,'fontsize', 12)
  xlabel('Time (s)','FontSize',12)
  ylabel('Q-axis current (A)','FontSize',12)
  legend('Measured','Estimated','Location','Northeast')
  grid on
  xlim([-0.1 0.1])
  ylim([-1 0])
  
 
  subplot(3,1,2)
  plot(t, sV(2,:), '-',t_Ts, xV(2,:), '--');
  set(gca,'fontsize', 12)
  xlabel('Time (s)','FontSize',12)
  ylabel('D-axis current (A)','FontSize',12)
  legend('Measured','Estimated','Location','Northeast')
  grid on
  xlim([-0.1 0.1])
  ylim([0 1])
  
  subplot(3,1,3)  
  plot(t, sV(3,:), '-',t_Ts, xV(3,:), '--');
  set(gca,'fontsize', 12)
  xlabel('Time (s)','FontSize',12)
  ylabel('Field current (A)','FontSize',12)
  legend('Measured','Estimated','Location','Southeast')
  grid on
  xlim([-0.1 0.1])
  ylim([0 .3])
  set(gcf,'Position',[0 0 850 1000])
  saveas(gcf,[nome 'Currents' '.png'])
  saveas(gcf,[nome 'Currents' '.fig'])
   
  figure(2)
  subplot(2,1,1)
  plot(t, sV(4,:), '-',t_Ts, xV(4,:), '--');
  set(gca,'fontsize', 12)
  xlabel('Time (s)','FontSize',12)
  ylabel('Angular velocity (rad/s)','FontSize',12)
  ylim([0 200])
  legend('Measured','Estimated','Location','Southeast')
  grid on
 
  
  subplot(2,1,2)
  plot(t, sV(5,:), '-',t_Ts, xV(5,:), '--');
  set(gca,'fontsize', 12)
  xlabel('Time (s)','FontSize',12)
  ylabel('Phase angle (rad)','FontSize',12)
  legend('Measured','Estimated','Location','Southeast')
  grid on
  set(gcf,'Position',[0 0 850 800])
  saveas(gcf,[nome 'Velocity-AngularPosition' '.png'])
  saveas(gcf,[nome 'Velocity-AngularPosition' '.fig'])

%   figure(4)
%   set(gcf,'Position',[0 0 850 800])
%   subplot(2,1,1)
%   plot(t_Ts, zV(1,:), 'k-',t_Ts, xV(1,:), 'g--');
%   set(gca,'fontsize', 12)
%   xlabel('Time (s)','FontSize',12)
%   ylabel('Q-axis current (A)','FontSize',12)
%   legend('Measured','Estimated')
%   grid on
%   
%   subplot(2,1,2)
%   plot(t_Ts, zV(2,:), 'k-',t_Ts, xV(2,:), 'g--');
%   set(gca,'fontsize', 12)
%   xlabel('Time (s)','FontSize',12)
%   ylabel('D-axis current (A)','FontSize',12)
%   legend('Measured','Estimated')
%   grid on
%   set(gcf,'Position',[0 0 850 800])
  
  nl = size(P_cov,1);
  na = size(P_cov,3);
  Pcov = zeros(nl,na);
  for k1 = 1:na
      for k2 = 1:nl
          Pcov(k2,k1) = P_cov(k2,k2,k1);
      end
  end

%   figure(5)
%   set(gcf,'Position',[0 0 850 340]);
%   plot(t_Ts,Pcov(1,:),'k-',t_Ts,Pcov(2,:),'k--',t_Ts,Pcov(3,:),'b-',t_Ts,Pcov(4,:),'b--',t_Ts,Pcov(5,:));
%   set(gca,'fontsize', 12)
%   xlabel('Time (s)','FontSize',12)
%   ylabel('State error covariance','FontSize',12)
%   legend('I_{q}','I_{d}','I_{fd}','\omega_m','\theta_m')
%   grid on
  
%    figure(5)
%    plot(t,0.1*ones(length(t),1),t_Ts,xV(8,:));%;,t_Ts,Smoothed_par);
%    xlabel('Time (s)')
%    ylabel('R_{s}')
%    legend('Actual','Estimated','Smoothed')
%    grid on
% 
%   figure(6)
%   plot(t,0.016*ones(length(t),1),t_Ts,xV(9,:));%;,t_Ts,Smoothed_par);
%   xlabel('Time (s)')
%   ylabel('R_{fd}')
%   legend('Actual','Estimated','Smoothed')
%   grid on
%  

if size(xV,1)>5
  figure(3)
  subplot(2,1,1)
%   subplot(2,1,1)
  plot(t,(0.680/0.85)*ones(length(t),1),t_Ts,xV(6,:)/0.85);%;,t_Ts,Smoothed_par);
  set(gca,'fontsize', 12)
  %set(gcf,'Position',[0 0 850 400])
  xlabel('Time (s)','FontSize',12)
  ylabel('L_{d}','FontSize',12)
  legend('Measured','Estimated','Location','Southeast')
  grid on
  ylim([0 1])
end
  
  if size(xV,1)>6
    subplot(2,1,2)
%   subplot(2,1,1)
  plot(t,0.571/0.85*ones(length(t),1),t_Ts,xV(7,:)/0.85);%;,t_Ts,Smoothed_par);
  set(gca,'fontsize', 12)
  %set(gcf,'Position',[0 0 850 400])
  xlabel('Time (s)','FontSize',12)
  ylabel('L_{q}','FontSize',12)
  legend('Measured','Estimated','Location','Southeast')
  grid on
  ylim([0 1])
  set(gcf,'Position',[0 0 850 1000])
  saveas(gcf,[nome 'Parameters' '.png'])
  saveas(gcf,[nome 'Parameters' '.fig'])
  end
  
%     if size(xV,1)>6
%     subplot(3,1,3)
%   plot(t,0.15*0.680*ones(length(t),1),t_Ts,0.15*xV(6,:));%;,t_Ts,Smoothed_par);
%   set(gca,'fontsize', 12)
%   %set(gcf,'Position',[0 0 850 400])
%   xlabel('Time (s)','FontSize',12)
%   ylabel('L_{ls}','FontSize',12)
%   legend('Measured','Estimated')
%   ylim([0 1])
%   grid on

% end
%   subplot(2,1,2)
%   plot(t,0.1*ones(length(t),1),t_Ts,xV(9,:));%;,t_Ts,Smoothed_par);
%   set(gca,'fontsize', 12)
%   %set(gcf,'Position',[0 0 850 400])
%   xlabel('Time (s)','FontSize',12)
%   ylabel('L_{md}','FontSize',12)
%   legend('Actual','Estimated','Smoothed')
%   grid on
%   set(gcf,'Position',[0 0 850 800])
  
%   figure(7)
%   set(gcf,'Position',[0 0 850 340])
%   plot(t_Ts,Pcov(8,:),t_Ts,Pcov(9,:));
%   set(gca,'fontsize', 12)
%   xlabel('Time (s)','FontSize',12)
%   ylabel('Covariance','FontSize',12)
%   legend('L_{mq}','L_{md}')
%   grid on
end
  