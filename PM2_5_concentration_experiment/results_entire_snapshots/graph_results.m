clear all, close all, clc;
%%
qiu = load('error_qui_batch.mat');
sob = load('error_sobolev_batch.mat');
nni = load('error_NNI.mat');
tik = load('error_tikhonov.mat');
%%
line_width = 1.5;
marker_size = 6;
font_size = 21;
width = 680;
heigth = 460;
path_figures = 'figures/';
mkdir(path_figures);
%%
semilogy(qiu.m,mean(qiu.RMSE),'ro--','LineWidth',line_width,'MarkerSize',marker_size);
hold on;
semilogy(nni.m,mean(nni.RMSE),'k*-','LineWidth',line_width,'MarkerSize',marker_size);
semilogy(tik.m,mean(tik.RMSE),'g-.','LineWidth',line_width,'MarkerSize',marker_size);
semilogy(sob.m,mean(sob.RMSE),'bs-','LineWidth',line_width,'MarkerSize',marker_size);
ylabel('RMSE','Interpreter','Latex');
xlabel('Snapshots sampling density','Interpreter','Latex');
xlim([sob.m(1) sob.m(end)]);
%ylim([2.9 6.05]);
lgd = legend({'Qiu','NNI','Tikhonov','GraphTRSS'},'Location','northeast');
lgd.NumColumns = 2;
set(lgd,'Interpreter','latex');
set(lgd,'color','none');
set(lgd,'Box','off');
title('b) Entire Snapshots Sampling','Interpreter','Latex');
get(gca);
set(gca,'FontName','times','FontSize',font_size,'TickLabelInterpreter','Latex');
set(gcf,'Position',[100,100,width,heigth]);
grid on;
saveas(gcf,[path_figures 'PM2_5_experiment.svg']);