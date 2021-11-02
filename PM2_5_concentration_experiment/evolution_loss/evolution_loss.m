clear all, close all, clc;
%%
addpath('../../functions');
%%
load('../graph_construction/full_graph.mat');
load('../PM2_5_concentration.mat');
x_matrix = myDataPM;
x_matrix = x_matrix(:,1:220);
%%
good_data = (x_matrix>0); % indicating the valid data
%%
signals_t = size(x_matrix,2);
%% load results
sob = load('../results/error_sobolev_batch.mat');
qiu = load('../results/error_qui_batch.mat');
%%
m = 0.1;  %Sampling density
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
        random_pattern = zeros(G.N,signals_t);
        for j=1:signals_t
            index_good_measures = find(good_data(:,j));
            sampled_index = index_good_measures(randperm(length(index_good_measures),...
                num_samples));
            random_pattern(sampled_index,j) = 1;
        end
        SampleMatrix = random_pattern;
        %% Random sampling
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