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
semilogy(anis.m(1:5),(anis.RMSE(1:5)),'ro--','LineWidth',line_width,'MarkerSize',marker_size);
hold on;
semilogy(chen.m(1:5),(chen.RMSE(1:5)),'k*-','LineWidth',line_width,'MarkerSize',marker_size);
semilogy(tsitsvero.m(1:5),(tsitsvero.RMSE(1:5)),'g-.','LineWidth',line_width,'MarkerSize',marker_size);
semilogy(random.m(1:5),mean(random.RMSE(:,1:5)),'bs-','LineWidth',line_width,'MarkerSize',marker_size);
ylabel('RMSE','Interpreter','Latex');
xlabel('Sampling density','Interpreter','Latex');
xlim([random.m(1) random.m(end-1)]);
lgd = legend({'Anis','Chen','Tsitsvero','Random'},'Location','southwest');
lgd.NumColumns = 2;
set(lgd,'Interpreter','latex');
set(lgd,'color','none');
set(lgd,'Box','off');
title('a) COVID-19 Global New Cases','Interpreter','Latex');
get(gca);
set(gca,'FontName','times','FontSize',font_size,'TickLabelInterpreter','Latex');
set(gcf,'Position',[100,100,width,heigth]);
grid on;
saveas(gcf,[path_figures 'covid_19_experiment.svg']);