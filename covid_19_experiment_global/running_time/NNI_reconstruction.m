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
%% Experiments
repetitions = 20;
running_time_NNI_random = zeros(repetitions,length(m));
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
        J = SampleMatrix;
        Y = J.*(x_matrix);
        tic
        x_recon = solver_NNI(J, Y, Position);
        running_time_NNI_random(ii,i) = toc;
    end
end
%%
results_path = '../results_time/';
mkdir(results_path);
save([results_path 'time_NNI_random.mat'],'running_time_NNI_random');