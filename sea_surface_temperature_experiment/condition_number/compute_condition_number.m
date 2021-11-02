clear all, close all, clc;
%%
load('../graph_construction/full_graph.mat');
load('../sea_surface_temperature.mat');
x_matrix = Data;
x_matrix = x_matrix(:,1:600);
%%
m = [0.5:0.1:0.9, 0.995];  %Sampling density
signals_t = size(x_matrix,2);
%%
load('../results/error_sobolev_batch.mat');
%%
lambda = alpha_set(best_alpha);
%epsilon = linspace(0.1,2);
epsilon = logspace(-4,7,70);
%beta = linspace(1e-3,6);
beta = [1];
%%
M = size(x_matrix,2);
Dh = temporal_difference_operator(M);
vec_Q = rand(M*G.N,1) > 0.5;
Q = spdiags(vec_Q,0,M*G.N,M*G.N);
%Q = sparse(diag(vec_Q));
condition_number_sob = zeros(length(epsilon),length(beta));
for i=1:length(epsilon)
    i
    for j=1:length(beta)
        condition_number_sob(i,j) = condest(Q+kron(lambda*Dh*Dh',(G.L+epsilon(i)*speye(G.N))^beta(j)));
    end
end
condition_number_lap = condest(Q+kron(lambda*Dh*Dh',G.L));
save(['condition_number_sob.mat'],'condition_number_sob','condition_number_lap','epsilon','beta');
semilogx(epsilon,condition_number_sob(:,1));
hold on;
size_epsilon = length(epsilon);
semilogx(epsilon,condition_number_lap*ones(size_epsilon,1));