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
%% load results
load('../results/sobolev_batch_parameter.mat');
%%
sobolev_matrix = G.L+epsilon_set(best_epsilon)*eye(G.N);
%% Symmetrization
sobolev_matrix = 0.5*(sobolev_matrix+sobolev_matrix');
sobolev_matrix = sparse(sobolev_matrix);
%%
param.niter = 40000;
param.beta = 0;
param.gamma = 0;
%% Experiment best parameters
repetitions = 100;
repetitions_sobolev = zeros(repetitions,length(m));
repetitions_qui = zeros(repetitions,length(m));
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
        param.J = SampleMatrix; % sampling matrix
        param.y = param.J.*(x_matrix); % sampled data
        param.T = TV_Temp(); % temporal difference operator 
        param.x0 = 0 * param.y;
        %% Sobolev
        param.L = sobolev_matrix; % Laplacian
        disp('Sobolev');
        [x_recon, iterations] = solver_BR_TVGS(param);
        repetitions_sobolev(ii,i) = iterations;
        %% Qui
        param.L = G.L; % Laplacian
        disp('Qiu');
        [x_recon, iterations] = solver_BR_TVGS(param);
        repetitions_qui(ii,i) = iterations;
    end
end
%%
results_path = '../results_convergence/';
mkdir(results_path);
save([results_path 'error_batch_random.mat'],'repetitions_sobolev','repetitions_qui');
%%
save([results_path 'm.mat'],'m');