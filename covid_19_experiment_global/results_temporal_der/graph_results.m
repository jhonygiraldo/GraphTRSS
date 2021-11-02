clear all, close all, clc;
%%
load('../results/error_sobolev_batch.mat')
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
semilogy(m,mean(RMSE),'ro--','LineWidth',line_width,'MarkerSize',marker_size);
hold on;
load('error_sobolev_batch_random_2.mat');
semilogy(m,mean(sqrt(error_matrix_sobolev_batch_random)),'k*-','LineWidth',line_width,'MarkerSize',marker_size);
load('error_sobolev_batch_random_3.mat');
semilogy(m,mean(sqrt(error_matrix_sobolev_batch_random)),'g-.','LineWidth',line_width,'MarkerSize',marker_size);
ylabel('RMSE','Interpreter','Latex');
xlabel('Sampling density','Interpreter','Latex');
xlim([m(1) m(end)]);
lgd = legend({'Step=1','Step=2','Step=3'},...
    'Location','northeast');
lgd.NumColumns = 2;
set(lgd,'Interpreter','latex');
set(lgd,'color','none');
set(lgd,'Box','off');
title('f) Temporal Difference Operators','Interpreter','Latex');
get(gca);
set(gca,'FontName','times','FontSize',font_size,'TickLabelInterpreter','Latex');
set(gcf,'Position',[100,100,width,heigth]);
grid on;
saveas(gcf,[path_figures 'covid_19_experiment.svg']);