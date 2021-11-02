clear all, close all, clc;
methods = {'NNI_random';'graph_reg';'tikhonov';'qui_batch';'puy';'sobolev_batch'};
for i=1:size(methods,1)
    load(['error_',methods{i},'.mat']);
    mean_RMSE = round(mean(mean(RMSE)),2);
    mean_MAE = round(mean(mean(MAE)),2);
    mean_MAPE = round(mean(mean(MAPE)),2);
    disp([methods{i},': $',num2str(mean_RMSE),'$ & $',num2str(mean_MAE),'$ & $',num2str(mean_MAPE),'$']);
end