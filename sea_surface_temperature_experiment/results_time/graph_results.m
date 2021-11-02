clear all, close all, clc;
%%
load('time_qui_batch_random.mat');
load('time_sobolev_batch.mat');
load('time_NNI_random.mat');
load('time_graph_reg_random.mat');
load('time_tikhonov_random.mat');
load('time_puy.mat');
load('m.mat');
%%
line_width = 1.5;
marker_size = 6;
font_size = 21;
width = 680;
heigth = 460;
path_figures = 'figures/';
mkdir(path_figures);
%%
semilogy(m,mean((running_time_matrix_qui_batch_random)),'ro--','LineWidth',line_width,'MarkerSize',marker_size);
hold on;
semilogy(m,mean((running_time_NNI_random)),'k*-','LineWidth',line_width,'MarkerSize',marker_size);
semilogy(m,mean((running_time_graph_reg_random)),'LineWidth',line_width,'MarkerSize',marker_size);
semilogy(m,mean((running_time_tikhonov_random)),'g-.','LineWidth',line_width,'MarkerSize',marker_size);
semilogy(m,mean((running_time_puy)),'m+:','LineWidth',line_width,'MarkerSize',marker_size);
semilogy(m,mean((running_time_sobolev_batch_random)),'bs-','LineWidth',line_width,'MarkerSize',marker_size);
ylabel('Time (s)','Interpreter','Latex');
xlabel('Sampling density','Interpreter','Latex');
xlim([m(1) m(end)]);
ylim([0.07 18])
lgd = legend({'Qiu','NNI','GR','Tikhonov','Puy','GraphTRSS'},'Location','northeast');
lgd.NumColumns = 3;
set(lgd,'Interpreter','latex');
set(lgd,'color','none');
set(lgd,'Box','off');
title('h) Running Time','Interpreter','Latex');
get(gca);
set(gca,'FontName','times','FontSize',font_size,'TickLabelInterpreter','Latex');
set(gcf,'Position',[100,100,width,heigth]);
grid on;
saveas(gcf,[path_figures 'sea_surface_experiment.svg']);