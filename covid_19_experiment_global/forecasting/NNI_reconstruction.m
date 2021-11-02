clear all, close all, clc;
%%
addpath('../../functions'); % Path to required functions
%%
load('../graph_construction/full_graph.mat');
load('../covid_19_new_cases.mat');
x_matrix = Data;
%%
predicted_snapshots = [1:10];
signals_t = size(x_matrix,2);
%% Experiments
repetitions = 1;
RMSE = zeros(repetitions,length(predicted_snapshots));
MAE = zeros(repetitions,length(predicted_snapshots));
MAPE = zeros(repetitions,length(predicted_snapshots));
for ii=1:repetitions
    ii
    for i=1:length(predicted_snapshots)
        %% Forecasting sampling
        forecasting_pattern = zeros(G.N,signals_t);
        forecasting_pattern(:,1:end-predicted_snapshots(i)) = 1;
        SampleMatrix = forecasting_pattern;
        %%
        J = SampleMatrix;
        Y = J.*(x_matrix);
        x_recon = solver_NNI(J, Y, Position);
        %%
        indx_non_sampled = find(SampleMatrix(:) == 0);
        x_vector_original = x_matrix(indx_non_sampled);
        x_vector_reconstructed = x_recon(indx_non_sampled);
        RMSE(ii,i) = sqrt(mean((x_vector_original-x_vector_reconstructed).^2));
        MAE(ii,i) = mean(abs(x_vector_original-x_vector_reconstructed));
        index_zero = find(x_vector_original == 0);
        x_vector_original(index_zero) = [];
        x_vector_reconstructed(index_zero) = [];
        MAPE(ii,i) = mean(abs((x_vector_original-x_vector_reconstructed)./x_vector_original));
    end
end
%%
results_path = '../results_forecasting/';
mkdir(results_path);
save([results_path 'error_NNI_random.mat'],'RMSE','MAE','MAPE','predicted_snapshots');