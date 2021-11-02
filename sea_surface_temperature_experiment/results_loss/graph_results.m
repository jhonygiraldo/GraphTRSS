clear all, close all, clc;
%%
load('loss_sob_qiu.mat')
%%
line_width = 1.5;
marker_size = 6;
font_size = 21;
width = 680;
heigth = 460;
path_figures = 'figures_unnormalized_Lap/';
mkdir(path_figures);
%%
max_iter = 50;
%% iterations
figure()
semilogy(mean(loss_qiu(:,1:max_iter)),'ro-','LineWidth',line_width,'MarkerSize',marker_size);
hold on;
semilogy(mean(loss_sobolev(:,1:max_iter)),'bs-','LineWidth',line_width,'MarkerSize',marker_size);
ylabel('Loss','Interpreter','Latex');
xlabel('Iteration','Interpreter','Latex');
lgd = legend({'Qiu','GraphTRSS'},'Location','northeast');
set(lgd,'Interpreter','latex');
set(lgd,'color','none');
set(lgd,'Box','off');
%title('g) Loss vs. Iterations','Interpreter','Latex');
title('c) Unnormalized Laplacian','Interpreter','Latex');
axis tight
grid on;
get(gca);
set(gca,'FontName','times','FontSize',font_size,'TickLabelInterpreter','Latex');
set(gcf,'Position',[100,100,width,heigth]);
saveas(gcf,[path_figures 'loss.svg']);