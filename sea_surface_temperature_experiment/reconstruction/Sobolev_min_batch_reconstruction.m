clear all, close all, clc;
%%
addpath('../../functions'); % Path to required functions
%%
load('../graph_construction/full_graph.mat');
load('../sea_surface_temperature.mat');
x_matrix = Data;
x_matrix = x_matrix(:,1:600);
%%
m = [0.1:0.1:0.9];  %Sampling density
signals_t = size(x_matrix,2);
%% alpha parameter
alpha_set = [1e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 1e2, 2e2, 5e2];
epsilon_set = [1e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 1e2, 2e2, 5e2];
error_sobolev_batch_random = cell(length(epsilon_set),length(alpha_set));
%%
sobolev_term_matrices = cell(1,length(epsilon_set));
for i=1:length(epsilon_set)
    laplacian = G.L+epsilon_set(i)*eye(G.N);
    %% Symmetrization
    sobolev_term = 0.5*(laplacian+laplacian.');
    sobolev_term_matrices{1,i} = sparse(sobolev_term);
end
%%
param.niter = 2000;
param.beta = 0;
param.gamma = 0;
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
        for jj=1:length(epsilon_set)
            param.L = sobolev_term_matrices{1,jj}; % Sobolev
            for j=1:length(alpha_set)
                param.alpha = alpha_set(j);
                param.x0 = 0 * param.y;
                x_recon = solver_BR_TVGS(param);
                indx_non_sampled = find(SampleMatrix(:) == 0);
                x_vector_original = x_matrix(indx_non_sampled);
                x_vector_reconstructed = x_recon(indx_non_sampled);
                error_sobolev_batch_random{jj,j}(ii,i) = mean((x_vector_original-x_vector_reconstructed).^2);
            end
        end
    end
end
%% Best parameters
mean_errors_sobolev_batch = zeros(size(error_sobolev_batch_random));
for i=1:size(error_sobolev_batch_random,1)
    for j=1:size(error_sobolev_batch_random,2)
        mean_errors_sobolev_batch(i,j) = mean(mean(error_sobolev_batch_random{i,j}));
    end
end
[best_epsilon, best_alpha] = find(mean_errors_sobolev_batch == min(min(mean_errors_sobolev_batch)));
%% Experiment best parameters
repetitions = 100;
RMSE = zeros(repetitions,length(m));
MAE = zeros(repetitions,length(m));
MAPE = zeros(repetitions,length(m));
param.niter = 20000;
param.alpha = alpha_set(best_alpha);
param.L = sobolev_term_matrices{1,best_epsilon}; % Sobolev
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
save([results_path 'error_sobolev_batch.mat'],'error_sobolev_batch_random',...
    'RMSE','MAE','MAPE','m','alpha_set','epsilon_set','best_epsilon','best_alpha');