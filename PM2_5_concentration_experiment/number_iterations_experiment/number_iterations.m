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