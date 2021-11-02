clear all, close all, clc;
methods = {'NNI_random';'graph_reg';'tikhonov';'qui_batch';'sobolev_batch'};
for i=1:size(methods,1)
    load(['error_',methods{i},'.mat']);
    best_RMSE = round(min(mean(RMSE)),2);
    best_MAE = round(min(mean(MAE)),2);
    best_MAPE = round(min(mean(MAPE)),2);
    disp([methods{i},': $',num2str(best_RMSE),'$ & $',num2str(best_MAE),'$ & $',num2str(best_MAPE),'$']);
end