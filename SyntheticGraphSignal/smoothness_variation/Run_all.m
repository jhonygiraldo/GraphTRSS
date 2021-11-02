clear; clc; close all;

%% run the methods
% natural neighbor interpolation
if ~exist('Error_NNI_RMSE.mat','file')
    Run_NNI
end

% graph regularization
if ~exist('Error_gsmooth_RMSE.mat','file')
    Run_graph_regularization
end

% graph-time Tikhonov
if ~exist('Error_station_RMSE.mat','file')
    Run_graph_time_Tikhonov
end

% Qiu batch method
if ~exist('Error_batch_RMSE.mat','file')
    Run_proposed_batch
end

% Puy method
if ~exist('Error_puy_RMSE.mat','file')
    Run_puy
end

% proposed Sobolev method
if ~exist('Error_Sobolev_RMSE.mat','file')
    Run_proposed_Sobolev
end

%% plot
load Error_NNI_RMSE
Error_NNI_RMSE = mean(Error_RMSE);

load Error_gsmooth_RMSE
Error_gsmooth_RMSE = mean(Error_RMSE);

load Error_station_RMSE
Error_station_RMSE = mean(Error_RMSE);

load Error_batch_RMSE
Error_batch_RMSE = mean(Error_RMSE);

load Error_puy_RMSE
Error_puy_RMSE = mean(Error_RMSE);

load Error_Sobolev_RMSE
Error_Sobolev_RMSE = mean(Error_RMSE);

line_width = 1.5;
marker_size = 6;
font_size = 21;
width = 680;
heigth = 460;
path_figures = 'figures/';
mkdir(path_figures);

figure; loglog(epsilon_set,Error_batch_RMSE,'ro--','LineWidth',line_width,'MarkerSize',marker_size);
hold on;
loglog(epsilon_set,Error_NNI_RMSE,'k*-','LineWidth',line_width,'MarkerSize',marker_size);
loglog(epsilon_set,Error_station_RMSE,'g-.','LineWidth',line_width,'MarkerSize',marker_size);
loglog(epsilon_set,Error_gsmooth_RMSE,'LineWidth',line_width,'MarkerSize',marker_size);
loglog(epsilon_set,Error_puy_RMSE,'m+:','LineWidth',line_width,'MarkerSize',marker_size);
loglog(epsilon_set,Error_Sobolev_RMSE,'bs-','LineWidth',line_width,'MarkerSize',marker_size);
ylabel('RMSE','Interpreter','Latex');
xlabel('Smoothness Level $\alpha$','Interpreter','Latex');
lgd = legend({'Qiu','NNI','Tikhonov','GR','Puy','GraphTRSS'},'Location','northwest');
lgd.NumColumns = 3;
set(lgd,'Interpreter','latex');
set(lgd,'color','none');
set(lgd,'Box','off');
title('c) Smoothness Variation','Interpreter','Latex');
get(gca);
set(gca,'FontName','times','FontSize',font_size,'TickLabelInterpreter','Latex');
set(gcf,'Position',[100,100,width,heigth]);
grid on;
axis tight
ylim([5e-2 120]);
saveas(gcf,[path_figures 'smoothness.svg']);