clear all, close all, clc;
%%
addpath('../../functions'); % Path to required functions
%%
load('../graph_construction/full_graph.mat');
load('../PM2_5_concentration.mat');
x_matrix = myDataPM;
x_matrix = x_matrix(:,1:220);
%%
good_data = (x_matrix>0); % indicating the valid data
%%
m = [0.1:0.05:0.45];  %Sampling density
signals_t = size(x_matrix,2);
%% alpha parameter
alpha_set = [1e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 1e2, 2e2, 5e2];
epsilon_set = [1e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 1e2, 2e2, 5e2];
beta_set = [0.5, 1.5, 2];
error_sobolev_batch_random = cell(length(epsilon_set),length(alpha_set),length(beta_set));
error_matrix_cell_sobolev_batch_random = cell(1,length(beta_set));
%%
load('../results/sobolev_batch_parameter.mat');
%%
sobolev_term_matrices = cell(length(beta_set),1);
for h=1:length(beta_set)
    laplacian = (G.L+epsilon_set(best_epsilon)*eye(G.N))^beta_set(h);
    %% Symmetrization
    sobolev_term = 0.5*(laplacian+laplacian.');
    sobolev_term_matrices{h,1} = sparse(sobolev_term);
end
%%
param.niter = 20000;
param.beta = 0;
param.gamma = 0;
%%
for h=1:length(beta_set)
    h
    %% Experiment best parameters
    repetitions = 10;
    error_matrix_sobolev_batch_random = zeros(repetitions,length(m));
    param.alpha = alpha_set(best_alpha);
    param.L = sobolev_term_matrices{h,1}; % Sobolev
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
            param.x0 = 0 * param.y;
            x_recon = solver_BR_TVGS(param);
            %%
            indx_non_sampled = find(SampleMatrix(:) == 0 & good_data(:) == 1);
            x_vector_original = x_matrix(indx_non_sampled);
            x_vector_reconstructed = x_recon(indx_non_sampled);
            error_matrix_sobolev_batch_random(ii,i) = mean(abs(x_vector_original-x_vector_reconstructed).^2);
        end
    end
    error_matrix_cell_sobolev_batch_random{h} = error_matrix_sobolev_batch_random;
    results_path = '../results_beta_/';
    mkdir(results_path);
    save([results_path 'error_sobolev_batch_random.mat'],'error_sobolev_batch_random','error_matrix_cell_sobolev_batch_random');
    %%
    save([results_path 'm.mat'],'m');
    save([results_path 'sobolev_batch_parameter.mat'],'alpha_set','epsilon_set','beta_set');
end
%%
results_path = '../results_beta/';
mkdir(results_path);
save([results_path 'error_sobolev_batch_random.mat'],'error_sobolev_batch_random','error_matrix_cell_sobolev_batch_random');
%%
save([results_path 'm.mat'],'m');
save([results_path 'sobolev_batch_parameter.mat'],'alpha_set','epsilon_set','beta_set');