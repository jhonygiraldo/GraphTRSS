clear all, close all, clc;
%%
qiu = load('error_qui_batch.mat');
sob = load('error_sobolev_batch.mat');
nni = load('error_NNI_random.mat');
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
semilogy(qiu.predicted_snapshots,(qiu.RMSE),'ro--','LineWidth',line_width,'MarkerSize',marker_size);
hold on;
semilogy(nni.predicted_snapshots,(nni.RMSE),'k*-','LineWidth',line_width,'MarkerSize',marker_size);
semilogy(tik.predicted_snapshots,(tik.RMSE),'g-.','LineWidth',line_width,'MarkerSize',marker_size);
semilogy(sob.predicted_snapshots,(sob.RMSE),'bs-','LineWidth',line_width,'MarkerSize',marker_size);
ylabel('RMSE','Interpreter','Latex');
xlabel('Temporal snapshots','Interpreter','Latex');
xlim([sob.predicted_snapshots(1) sob.predicted_snapshots(end)]);
lgd = legend({'Qiu','NNI','Tikhonov','GraphTRSS'},'Location','northwest');
lgd.NumColumns = 2;
set(lgd,'Interpreter','latex');
set(lgd,'color','none');
set(lgd,'Box','off');
title('c) Forecasting','Interpreter','Latex');
get(gca);
set(gca,'FontName','times','FontSize',font_size,'TickLabelInterpreter','Latex');
set(gcf,'Position',[100,100,width,heigth]);
grid on
saveas(gcf,[path_figures 'covid_global_experiment.svg']);