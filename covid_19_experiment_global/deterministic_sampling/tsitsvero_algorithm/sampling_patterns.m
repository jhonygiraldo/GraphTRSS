clear all, close all, clc;
load('../../graph_construction/full_graph.mat');
G = gsp_compute_fourier_basis(G);
%%
N = G.N;
m = [0.5:0.1:0.9,0.995];  %Sampling density
m = round(N*m);
m_max = max(m);
sampling_index_tsitsvero = zeros(m_max);
bandwidth = 155;
for(i=1:length(m))
    i
    [~,sampling_index_tsitsvero] = sampling_tsitsvero_algorithm(G.U,m(i),bandwidth);
end
sampling_patterns_tsitsvero = zeros(length(m),N);
for i=1:length(m)
    sampling_patterns_tsitsvero(i,sampling_index_tsitsvero(1:m(i))) = 1;
end
save('sampling_patterns_tsitsvero.mat','sampling_patterns_tsitsvero');