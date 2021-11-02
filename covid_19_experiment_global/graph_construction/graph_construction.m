clear all, close all, clc;
%%
load('../covid_19_new_cases.mat');
points = Position;
%%
N = size(points,1);
knn_param = 10;
[Idx Dist] = knnsearch(points,points,'K',knn_param+1);
sigma = mean(mean(Dist));
W = spalloc(N,N,(2*N*knn_param));
for i=1:N
    i/N
    W(i,Idx(i,2:end)) = exp(-(Dist(i,2:end).^2)./(sigma^2));
    W(Idx(i,2:end),i) = W(i,Idx(i,2:end));
end
%%
G.N = N;
G.W = W;
G.coords = points;
G.type = 'nearest neighbors';
G.sigma = sigma;
G = gsp_graph_default_parameters(G);
G = gsp_estimate_lmax(G);
save(['full_graph.mat'],'G','Idx','Dist');
plot(graph(G.W),'XData',G.coords(:,2),'YData',G.coords(:,1));