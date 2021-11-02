clear all, close all, clc;
%%
load('error_sobolev_batch_random.mat')
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
semilogy(m,mean(sqrt(error_matrix_cell_sobolev_batch_random{1,1})),'ro--','LineWidth',line_width,'MarkerSize',marker_size);
hold on;
load('../results/error_sobolev_batch.mat');
semilogy(m,mean((RMSE)),'k*-','LineWidth',line_width,'MarkerSize',marker_size);
load('error_sobolev_batch_random_2.mat')
semilogy(m,mean(sqrt(error_matrix_cell_sobolev_batch_random{1,2})),'g-.','LineWidth',line_width,'MarkerSize',marker_size);
semilogy(m,mean(sqrt(error_matrix_cell_sobolev_batch_random{1,3})),'bs-','LineWidth',line_width,'MarkerSize',marker_size);
ylabel('RMSE','Interpreter','Latex');
xlabel('Sampling density','Interpreter','Latex');
xlim([m(1) m(end)]);
lgd = legend({'$\beta=0.5$','$\beta=1$','$\beta=1.5$','$\beta=2$'},...
    'Location','northeast');
lgd.NumColumns = 2;
set(lgd,'Interpreter','latex');
set(lgd,'color','none');
set(lgd,'Box','off');
title('e) Variation $\beta$','Interpreter','Latex');
get(gca);
set(gca,'FontName','times','FontSize',font_size,'TickLabelInterpreter','Latex');
set(gcf,'Position',[100,100,width,heigth]);
grid on;
saveas(gcf,[path_figures 'sea_surface_experiment.svg']);