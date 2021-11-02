clear all, close all, clc;
%%
addpath('../../functions'); % Path to required functions
%%
load('../graph_construction/full_graph.mat');
load('../PM2_5_concentration.mat');
x_matrix = myDataPM;
x_matrix = x_matrix(:,1:220);
%%
good_data = (x_matrix > 0); % indicating the valid data
%%
m = [0.1:0.05:0.45];  %Sampling density
signals_t = size(x_matrix,2);
%% alpha parameter
beta_set = [1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50];
error_graph_reg_random = cell(1,length(beta_set));
%%
param.L = G.L; % Graph Laplacian
param.niter = 2000;
param.alpha = 0;
param.gamma = 0;
%%
repetitions = 1;
for ii=1:repetitions
    ii
    for i=1:length(m)
        num_samples = round(m(i)*G.N);
        %% Random sampling
        random_pattern = zeros(G.N,signals_t);
        for j=1:signals_t
            index_good_measures = find(good_data(:,j));
            sampled_index = index_good_measures(randperm(length(index_good_measures),...
                num_samples));
            random_pattern(sampled_index,j) = 1;
        end
        SampleMatrix = random_pattern;
        %%
        param.J = SampleMatrix; % sampling matrix
        param.y = param.J.*(x_matrix); % sampled data
        param.T = TV_Temp(); % temporal difference operator
        for j=1:length(beta_set)
            param.beta = beta_set(j);
            param.x0 = 0 * param.y;
            x_recon = solver_BR_TVGS(param);
            %%
            indx_non_sampled = find(SampleMatrix(:) == 0 & good_data(:) == 1);
            x_vector_original = x_matrix(indx_non_sampled);
            x_vector_reconstructed = x_recon(indx_non_sampled);
            error_graph_reg_random{1,j}(ii,i) = mean((x_vector_original-x_vector_reconstructed).^2);
        end
    end
end
%% Best parameters
mean_errors_graph_reg = zeros(size(error_graph_reg_random));
for i=1:size(error_graph_reg_random,2)
    mean_errors_graph_reg(1,i) = mean(error_graph_reg_random{1,i});
end
[best_x, best_beta] = find(mean_errors_graph_reg == min(mean_errors_graph_reg));
%% Experiment best parameters
repetitions = 100;
RMSE = zeros(repetitions,length(m));
MAE = zeros(repetitions,length(m));
MAPE = zeros(repetitions,length(m));
param.niter = 20000;
%%
param.alpha = 0;
param.gamma = 0;
for ii=1:repetitions
    ii
    for i=1:length(m)
        num_samples = round(m(i)*G.N);
        %% Random sampling
        random_pattern = zeros(G.N,signals_t);
        for j=1:signals_t
            index_good_measures = find(good_data(:,j));
            sampled_index = index_good_measures(randperm(length(index_good_measures),...
                num_samples));
            random_pattern(sampled_index,j) = 1;
        end
        SampleMatrix = random_pattern;
        %%
        param.J = SampleMatrix; % sampling matrix
        param.y = param.J.*(x_matrix); % sampled data
        param.T = TV_Temp(); % temporal difference operator
        param.beta = beta_set(best_beta);
        param.x0 = 0 * param.y;
        x_recon = solver_BR_TVGS(param);
        %%
        indx_non_sampled = find(SampleMatrix(:) == 0 & good_data(:) ==1);
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
save([results_path 'error_graph_reg_random.mat'],'error_graph_reg_random',...
	'RMSE','MAE','MAPE','m','beta_set','best_beta');