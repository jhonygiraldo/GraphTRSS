clear all, close all, clc;
load('../../graph_construction/full_graph.mat');
%%
N = G.N;
m = [0.1:0.1:0.9];  %Sampling density
m = round(N*m);
kpow=7;
L=diag(G.d)-G.W;
L_kpow=L^kpow;
sampling_patterns_anis = zeros(length(m),N);
for(i=1:length(m))
    i
    [sampling_patterns_anis(i,:),~,~]=compute_opt_set_inc(L_kpow, kpow, m(i));
end
save('sampling_patterns_anis.mat','sampling_patterns_anis');