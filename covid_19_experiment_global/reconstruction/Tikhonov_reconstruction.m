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
%% alpha parameter
beta_set = [1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 1e2, 2e2, 5e2];
gamma_set = [1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1];
error_tikhonov_random = cell(length(beta_set),length(gamma_set));
%%
param.niter = 2000;
param.gamma = 0;
param.L = G.L; % Laplacian
%%
repetitions = 1;
for ii=1:repetitions
    ii
    for i=1:length(m)
        num_samples = round(m(i)*G.N);
        %% Random sampling
        random_pattern = zeros(signals_t,G.N);
        for j=1:signals_t
            random_pattern(j,randperm(G.N,num_samples)) = 1;
        end
        SampleMatrix = random_pattern';
        param.J = SampleMatrix; % sampling matrix
        param.y = param.J.*(x_matrix); % sampled data
        param.T = TV_Temp(); % temporal difference operator
        param.alpha = 0;
        for jj=1:length(beta_set)
            param.beta = beta_set(jj);
            for j=1:length(gamma_set)
                param.gamma = gamma_set(j);
                param.x0 = 0 * param.y;
                x_recon = solver_BR_TVGS(param);
                indx_non_sampled = find(SampleMatrix(:) == 0);
                x_vector_original = x_matrix(indx_non_sampled);
                x_vector_reconstructed = x_recon(indx_non_sampled);
                error_tikhonov_random{jj,j}(ii,i) = mean((x_vector_original-x_vector_reconstructed).^2);
            end
        end
    end
end
%% Best parameters
mean_errors_tikhonov = zeros(size(error_tikhonov_random));
for i=1:size(error_tikhonov_random,1)
    for j=1:size(error_tikhonov_random,2)
        mean_errors_tikhonov(i,j) = mean(mean(error_tikhonov_random{i,j}));
    end
end
[best_beta, best_gamma] = find(mean_errors_tikhonov == min(min(mean_errors_tikhonov)));
%% Experiment best parameters
repetitions = 100;
RMSE = zeros(repetitions,length(m));
MAE = zeros(repetitions,length(m));
MAPE = zeros(repetitions,length(m));
param.niter = 20000;
param.beta = beta_set(best_beta);
param.gamma = gamma_set(best_gamma);
param.alpha = 0;
param.L = G.L; % Laplacian
for ii=1:repetitions
    ii
    for i=1:length(m)
        num_samples = round(m(i)*G.N);
        %% Random sampling
        random_pattern = zeros(signals_t,G.N);
        for j=1:signals_t
            random_pattern(j,randperm(G.N,num_samples)) = 1;
        end
        SampleMatrix = random_pattern';
        param.J = SampleMatrix; % sampling matrix
        param.y = param.J.*(x_matrix); % sampled data
        param.T = TV_Temp(); % temporal difference operator
        param.x0 = 0 * param.y;
        x_recon = solver_BR_TVGS(param);
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
results_path = '../results/';
mkdir(results_path);
save([results_path 'error_tikhonov.mat'],'error_tikhonov_random',...
    'RMSE','MAE','MAPE','m','beta_set','gamma_set','best_beta','best_gamma');