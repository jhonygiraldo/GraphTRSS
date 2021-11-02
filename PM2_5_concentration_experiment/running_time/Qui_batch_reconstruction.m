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
%%
load('../results/qui_batch_parameter.mat');
%%
param.L = G.L; % Graph Laplacian
param.niter = 20000;
param.beta = 0;
param.gamma = 0;
%% Experiment best parameters
repetitions = 20;
running_time_matrix_qui_batch_random = zeros(repetitions,length(m));
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
        param.alpha = alpha_set(best_alpha);
        param.x0 = 0 * param.y;
        tic
        x_recon = solver_BR_TVGS(param);
        running_time_matrix_qui_batch_random(ii,i) = toc;
    end
end
%%
results_path = '../results_time/';
mkdir(results_path);
save([results_path 'time_qui_batch_random.mat'],'running_time_matrix_qui_batch_random');
%%
save([results_path 'm.mat'],'m');