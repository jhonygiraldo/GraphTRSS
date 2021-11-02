clear all, close all, clc;
%%
load('../PM2_5_concentration.mat');
points = Position;
%%
N = size(points,1);
knn_param = 10;
[Idx Dist] = knnsearch(points,points,'K',knn_param+1);
sigma = mean(mean(Dist));
W = sparse(N,N);
for i=1:N
    for j=2:knn_param+1
        W(i,Idx(i,j)) = exp(-(Dist(i,j)^2)/(sigma^2));
        W(Idx(i,j),i) = W(i,Idx(i,j));
    end
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