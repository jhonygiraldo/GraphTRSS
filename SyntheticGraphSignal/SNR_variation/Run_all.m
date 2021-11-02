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

% proposed Puy method
if ~exist('Error_puy_RMSE.mat','file')
    Run_Puy
end

% proposed Sobolev method
if ~exist('Error_Sobolev_RMSE.mat','file')
    Run_proposed_Sobolev
end

%% plot
load Error_NNI_RMSE
Error_NNI_RMSE = Error_RMSE;

load Error_gsmooth_RMSE
Error_gsmooth_RMSE = Error_RMSE;

load Error_station_RMSE
Error_station_RMSE = Error_RMSE;

load Error_batch_RMSE
Error_batch_RMSE = Error_RMSE;

load Error_puy_RMSE
Error_puy_RMSE = Error_RMSE;

load Error_Sobolev_RMSE
Error_Sobolev_RMSE = Error_RMSE;

load paramAWD;
[N,T] = size(Temp);
noise_set2 = -20*log10(noise_set* sqrt(N*T) ./ norm(Temp));
Error_NNI_RMSE2 = -20*log10(Error_NNI_RMSE * sqrt(N*T) ./ norm(Temp));
Error_gsmooth_RMSE2 = -20*log10(Error_gsmooth_RMSE * sqrt(N*T) ./ norm(Temp));
Error_station_RMSE2 = -20*log10(Error_station_RMSE * sqrt(N*T) ./ norm(Temp));
Error_batch_RMSE2 = -20*log10(Error_batch_RMSE * sqrt(N*T) ./ norm(Temp));
Error_puy_RMSE2 = -20*log10(Error_puy_RMSE * sqrt(N*T) ./ norm(Temp));
Error_Sobolev_RMSE2 = -20*log10(Error_Sobolev_RMSE * sqrt(N*T) ./ norm(Temp));

line_width = 1.5;
marker_size = 6;
font_size = 21;
width = 680;
heigth = 460;
path_figures = 'figures/';
mkdir(path_figures);

figure; semilogy(noise_set2,Error_batch_RMSE2,'ro--','LineWidth',line_width,'MarkerSize',marker_size);
hold on;
semilogy(noise_set2,Error_NNI_RMSE2,'k*-','LineWidth',line_width,'MarkerSize',marker_size);
semilogy(noise_set2,Error_station_RMSE2,'g-.','LineWidth',line_width,'MarkerSize',marker_size);
semilogy(noise_set2,Error_gsmooth_RMSE2,'LineWidth',line_width,'MarkerSize',marker_size);
semilogy(noise_set2,Error_puy_RMSE2,'m+:','LineWidth',line_width,'MarkerSize',marker_size);
semilogy(noise_set2,Error_Sobolev_RMSE2,'bs-','LineWidth',line_width,'MarkerSize',marker_size);
ylabel('Output SNR (dB)','Interpreter','Latex');
xlabel('Input SNR (dB)','Interpreter','Latex');
lgd = legend({'Qiu','NNI','Tikhonov','GR','Puy','GraphTRSS'},'Location','best');
lgd.NumColumns = 3;
set(lgd,'Interpreter','latex');
set(lgd,'color','none');
set(lgd,'Box','off');
title('b) SNR Variation','Interpreter','Latex');
get(gca);
set(gca,'FontName','times','FontSize',font_size,'TickLabelInterpreter','Latex');
set(gcf,'Position',[100,100,width,heigth]);
grid on;
axis tight
ylim([5 38]);
saveas(gcf,[path_figures 'SNR.svg']);