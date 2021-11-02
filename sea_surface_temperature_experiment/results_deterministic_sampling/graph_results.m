clear all, close all, clc;
%%
random = load('../results/error_sobolev_batch.mat');
anis = load('error_sobolev_batch_anis.mat');
chen = load('error_sobolev_batch_chen.mat');
tsitsvero = load('error_sobolev_batch_tsitsvero.mat');
%%
line_width = 1.5;
marker_size = 6;
font_size = 21;
width = 680;
heigth = 460;
path_figures = 'figures/';
mkdir(path_figures);
%%
semilogy(anis.m,(anis.RMSE),'ro--','LineWidth',line_width,'MarkerSize',marker_size);
hold on;
semilogy(chen.m,(chen.RMSE),'k*-','LineWidth',line_width,'MarkerSize',marker_size);
semilogy(tsitsvero.m,(tsitsvero.RMSE),'g-.','LineWidth',line_width,'MarkerSize',marker_size);
semilogy(random.m,mean(random.RMSE),'bs-','LineWidth',line_width,'MarkerSize',marker_size);
ylabel('RMSE','Interpreter','Latex');
xlabel('Sampling density','Interpreter','Latex');
xlim([random.m(1) random.m(end)]);
lgd = legend({'Anis','Chen','Tsitsvero','Random'},'Location','northwest');
lgd.NumColumns = 2;
set(lgd,'Interpreter','latex');
set(lgd,'color','none');
set(lgd,'Box','off');
title('b) Sea Surface Temperature','Interpreter','Latex');
get(gca);
set(gca,'FontName','times','FontSize',font_size,'TickLabelInterpreter','Latex');
set(gcf,'Position',[100,100,width,heigth]);
grid on;
saveas(gcf,[path_figures 'sea_experiment.svg']);