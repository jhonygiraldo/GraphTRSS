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
%% load results
load('../results/sobolev_batch_parameter.mat');
%%
epsilon_set = [1e-5, 1e-4, 1e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10];
%%
sobolev_term_matrices = cell(1,length(epsilon_set));
for i=1:length(epsilon_set)
    laplacian = G.L+epsilon_set(i)*eye(G.N);
    %% Symmetrization
    sobolev_term = 0.5*(laplacian+laplacian.');
    sobolev_term_matrices{1,i} = sparse(sobolev_term);
end
%%
param.niter = 40000;
param.beta = 0;
param.gamma = 0;
%% Experiment
repetitions = 10;
repetitions_sobolev = cell(repetitions,length(epsilon_set));
error_matrix_sobolev = cell(repetitions,length(epsilon_set));
param.alpha = alpha_set(best_alpha);
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
        %%
        param.J = SampleMatrix; % sampling matrix
        param.y = param.J.*(x_matrix); % sampled data
        param.T = TV_Temp(); % temporal difference operator 
        param.x0 = 0 * param.y;
        for j=1:length(epsilon_set)
            %% Sobolev
            param.L = sobolev_term_matrices{1,j}; % Laplacian
            [x_recon, iterations] = solver_BR_TVGS(param);
            repetitions_sobolev{ii,j}(1,i) = iterations;
            %%
            indx_non_sampled = find(SampleMatrix(:) == 0);
            x_vector_original = x_matrix(indx_non_sampled);
            x_vector_reconstructed = x_recon(indx_non_sampled);
            error_matrix_sobolev{ii,j}(1,i) = mean(abs(x_vector_original-x_vector_reconstructed).^2);
        end
    end
end
%%
results_path = '../results_convergence/';
mkdir(results_path);
save([results_path 'iterations_vs_epsilon.mat'],'repetitions_sobolev','error_matrix_sobolev');
%%
save([results_path 'epsilon_set.mat'],'epsilon_set');