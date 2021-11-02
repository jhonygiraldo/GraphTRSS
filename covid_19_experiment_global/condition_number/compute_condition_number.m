clear all, close all, clc;
%%
load('../graph_construction/full_graph.mat');
load('../covid_19_new_cases.mat');
x_matrix = Data;
%%
m = [0.5:0.1:0.9, 0.995];  %Sampling density
signals_t = size(x_matrix,2);
%%
load('../results/error_sobolev_batch.mat');
%%
lambda = alpha_set(best_alpha);
epsilon = logspace(-4,4);
beta = [1:6];
%%
M = size(x_matrix,2);
Dh = sparse(temporal_difference_operator(M));
vec_Q = rand(M*G.N,1) > 0.5;
Q = sparse(diag(vec_Q));
condition_number_sob = zeros(length(epsilon),length(beta));
for i=1:length(epsilon)
    i
    for j=1:length(beta)
        condition_number_sob(i,j) = condest(Q+kron(lambda*Dh*Dh',(G.L+epsilon(i)*speye(G.N))^beta(j)));
    end
end
condition_number_lap = condest(Q+kron(lambda*Dh*Dh',G.L));
save(['condition_number_sob.mat'],'condition_number_sob','condition_number_lap','epsilon','beta');