clear all, close all, clc;
%%
load('../paramAWD.mat');
x_matrix = Temp;
%%
lambda = 0.1;
epsilon = logspace(-1,1);
beta = [0.5, 1, 1.5];
%%
[N,M] = size(x_matrix);
Dh = temporal_difference_operator(M);
vec_Q = rand(M*N,1) > 0.5;
Q = spdiags(vec_Q,0,M*N,M*N);
condition_number_sob = zeros(length(epsilon),length(beta));
for i=1:length(epsilon)
    i
    for j=1:length(beta)
        condition_number_sob(i,j) = condest(Q+kron(lambda*Dh*Dh',(L+epsilon(i)*speye(N))^beta(j)));
    end
end
condition_number_lap = condest(Q+kron(lambda*Dh*Dh',L));
save(['condition_number_sob.mat'],'condition_number_sob','condition_number_lap','epsilon','beta');