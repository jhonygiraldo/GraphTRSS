clear all, close all, clc;
%%
addpath('../../functions'); % Path to required functions
%%
load('../graph_construction/full_graph.mat');
load('../covid_19_new_cases.mat');
x_matrix = Data;
%%
m = [0.5:0.1:0.9, 0.995];  %Sampling density
signals_t = size(x_matrix,2);
%% Experiments
repetitions = 100;
RMSE = zeros(repetitions,length(m));
MAE = zeros(repetitions,length(m));
MAPE = zeros(repetitions,length(m));
for ii=1:repetitions
    ii
    for i=1:length(m)
        num_snapshots = round(m(i)*signals_t);
        %% Random sampling
        random_pattern = zeros(G.N,signals_t);
        random_pattern(:,randperm(signals_t,num_snapshots)) = 1;
        SampleMatrix = random_pattern;
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
results_path = '../results_entire_snapshots/';
mkdir(results_path);
save([results_path 'error_NNI_random.mat'],'RMSE','MAE','MAPE','m');