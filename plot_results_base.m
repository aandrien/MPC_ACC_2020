%% Default code
clc; close all; clear;
saveplot = 0;

%% Main code
data_MPC = load('Sim_V1_MPC_step.mat');
data_base = load('Sim_V1_MPC_step_base.mat');
% load('Sim_V1_MPC.mat')

p_MPC = data_MPC.pos;
v_MPC = data_MPC.vel;
a_MPC = data_MPC.acc;

p_base = data_base.pos;
v_base = data_base.vel;
a_base = data_base.acc;

L = data_MPC.L;

t = data_MPC.tout;
tend = data_MPC.tend;
plot_t = 1:3:tend/data_MPC.dt;
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

lineColors = linspecer(2);
lineStyles = {'-','--','-.'};

figure(1)
subplot(3,1,1)
plot(t(plot_t),p_MPC((plot_t),1),'color',lineColors(1,:),'LineStyle',lineStyles{1})
hold on
plot(t(plot_t),p_base((plot_t),1),'color',lineColors(2,:),'LineStyle',lineStyles{2})

legend('MPC','Low-gain')
ylabel('Position $[m]$','interpreter','latex')

subplot(3,1,2)
plot(t(plot_t),p_MPC((plot_t),1),'color',lineColors(1,:),'LineStyle',lineStyles{1})
hold on
plot(t(plot_t),p_base((plot_t),1),'color',lineColors(2,:),'LineStyle',lineStyles{2})
ylabel('Position Error $[m]$','interpreter','latex')

subplot(3,1,3)
plot(t(plot_t),a_MPC((plot_t),1),'color',lineColors(1,:),'LineStyle',lineStyles{1})
hold on
plot(t(plot_t),a_base((plot_t),1),'color',lineColors(2,:),'LineStyle',lineStyles{2})
plot([t(1) t(plot_t(end))],[L L],'k--')
plot([t(1) t(plot_t(end))],[-L -L],'k--')
ylim([-1.1*L 1.1*L])
ylabel('Acceleration $[\frac{m}{s^2}]$','interpreter','latex')
xlabel('Time $[s]$','interpreter','latex')

if(saveplot)
    set(gcf,'color','w');
    export_fig ..\..\V12\Pictures\MPC_step_base.eps
end


