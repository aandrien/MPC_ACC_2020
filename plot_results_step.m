%% Default code
clc; close all; clear;
saveplot = 0;

%% Main code
load('Sim_V1_MPC_step.mat');

p_MPC = pos;
v_MPC = vel;
a_MPC = acc;

t = tout;
plot_t = 1:3:t(end)/dt;
%% Figures
% figure_configuration_IEEE_standard
set(0, 'DefaultLineLineWidth', 2);
set(0,'DefaultAxesXGrid','on')
set(0,'DefaultAxesYGrid','on')
set(0,'defaultAxesFontName','Times New Roman');
set(0,'defaultAxesFontSize',12);
set(0,'defaultTextFontName','Times New Roman');
set(0,'defaultTextFontSize',12);
set(0,'defaultLegendFontName','Times New Roman');
set(0,'defaultLegendFontSize',12);

lineColors = linspecer(4);
lineStyles = {'-','--','-.'};

figure(1)
subplot(3,1,1)
hold on
for ii = 1:3
    plot(t(plot_t),p_MPC((plot_t),ii),'color',lineColors(ii,:),'LineStyle',lineStyles{ii})
end
legend({'$x$','$y$','$z$'},'interpreter','latex')
ylabel('Position $[m]$','interpreter','latex')

subplot(3,1,2)
hold on
for ii = 1:3
    plot(t(plot_t),v_MPC((plot_t),ii),'color',lineColors(ii,:),'LineStyle',lineStyles{ii})
end
ylabel('Velocity $[\frac{m}{s}]$','interpreter','latex')

subplot(3,1,3)
hold on
for ii = 1:3
    plot(t(plot_t),a_MPC((plot_t),ii),'color',lineColors(ii,:),'LineStyle',lineStyles{ii})
end
plot([t(1) t(plot_t(end))],[L L],'k--')
plot([t(1) t(plot_t(end))],[-L -L],'k--')
ylim([-1.01*L 1.01*L])
ylabel('Acceleration $[\frac{m}{s^2}]$','interpreter','latex')
xlabel('Time $[s]$','interpreter','latex')

if(saveplot)
    set(gcf,'color','w');
    export_fig ..\..\V12\Pictures\MPC_step.eps
end

