clear all, close all, clc;
load('../../graph_construction/full_graph.mat');
G = gsp_compute_fourier_basis(G);
%%
N = G.N;
m = [0.1:0.1:0.9];  %Sampling density
m = round(N*m);
m_max = max(m);
K = 23;
V_K = G.U(:,1:K);
[sampling_operator, sampled_indexes] = sampling_chen_algorithm(V_K,m_max,N);
sampling_operator = zeros(length(m),N);
for i=1:length(m)
    sampling_operator(i,sampled_indexes(1:m(i))) = 1;
end
save('sampling_operator_chen.mat','sampling_operator');