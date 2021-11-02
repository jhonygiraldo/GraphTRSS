clear all, close all, clc;
%%
addpath('../../functions');
%%
load('../graph_construction/full_graph.mat');
load('../covid_19_new_cases.mat');
x_matrix = Data;
%%
signals_t = size(x_matrix,2);
%% load results
sob = load('../results/error_sobolev_batch.mat');
qiu = load('../results/error_qui_batch.mat');
%%
m = 0.5;  %Sampling density
%%
sobolev_matrix = G.L+sob.epsilon_set(sob.best_epsilon)*eye(G.N);
%% Symmetrization
sobolev_matrix = 0.5*(sobolev_matrix+sobolev_matrix');
sobolev_matrix = sparse(sobolev_matrix);
%%
param.niter = 10000;
param.beta = 0;
param.gamma = 0;
%% Experiment best parameters
repetitions = 50;
loss_sobolev = zeros(repetitions,param.niter);
loss_qiu = zeros(repetitions,param.niter);
repetitions_sobolev = zeros(repetitions,length(m));
repetitions_qui = zeros(repetitions,length(m));
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
        param.alpha = sob.alpha_set(sob.best_alpha);
        disp('Sobolev');
        [x_recon, iterations, loss_sobolev(ii,:)] = solver_BR_TVGS(param);
        repetitions_sobolev(ii,i) = iterations;
        %% Qui
        param.L = G.L; % Laplacian
        param.alpha = qiu.alpha_set(qiu.best_alpha);
        disp('Qiu');
        [x_recon, iterations, loss_qiu(ii,:)] = solver_BR_TVGS(param);
        repetitions_qui(ii,i) = iterations;
    end
end
%%
results_path = '../results_loss/';
mkdir(results_path);
save([results_path 'loss_sob_qiu.mat'],'loss_sobolev','loss_qiu',...
    'repetitions_sobolev','repetitions_qui','m');