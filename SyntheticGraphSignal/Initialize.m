% ========================================================================
% Time-varying Graph Signal Reconstruction,
%
% Copyright(c) 2017 Kai Qiu and Yuantao Gu
% All Rights Reserved.
% ----------------------------------------------------------------------
% 100 points are generated randomly.
% The dimension of the synthetic dataset is 100x600. 
% 
% Version 1.0
% Written by Kai Qiu (q1987k@163.com)
%----------------------------------------------------------------------

clear; clc; close all;
addpath('../../utilities');


%% graph construction
N = 100;
Position = zeros(N, 2);
point_select = randperm(N*N, N); % selected points
for i_point = 1 : N
    Position(i_point, 1) = mod(point_select(i_point)-1, N) +1; % x coordinate
    Position(i_point, 2) = floor((point_select(i_point)-1)/N) + 1; % y coordinate
end

Dist = zeros(N); % distances between each two points
for i = 1:N 
    for j = i+1:N
        Dist(i,j) = norm(Position(i,:) - Position(j,:));
        Dist(j,i) = Dist(i,j);
    end
end
A=zeros(N); % Adjacency matrix
W=zeros(N); % Weighted matrix
k=5; % k-NN method
for i=1:N
    [sortd,ind]=sort(Dist(i,:),'ascend');
    A(i,ind(2:k+1))=1;
    A(ind(2:k+1),i)=1;
    W(i,ind(2:k+1))=1./(Dist(i,ind(2:k+1))).^2;
    W(ind(2:k+1),i)=W(i,ind(2:k+1));
end
W=W/max(max(W));
D=diag(sum(W));
L=D-W;

%% plot the graph
C = linspecer(7);
figure; hold on;
plot(Position(:,1), Position(:,2),'o', 'color', C(4,:), 'Markerfacecolor', C(4,:),'Markersize', 4);
[ki,kj]=find(A);
plot([Position(ki,1)';Position(kj,1)'],[Position(ki,2)';Position(kj,2)'],'LineWidth',1, 'color', C(2,:));


%% generate the time-varying graph signal
[V,Lambda]=eig(L);
Lambda(1,1)=0;
lambda = diag(Lambda);
lambdaHalfInv = 1./sqrt(lambda);
lambdaHalfInv(1) = 0;
LHalfInv = V*diag(lambdaHalfInv)*V';

T = 600;
Temp = zeros(N, T);
ftmp = V'* randn(N,1);
ftmp(11:end) = ftmp(11:end)/100; % Weakened the high frequency components
ftmp = V*ftmp;
Temp(:,1) = ftmp / norm(ftmp) * 100;

epsilon = 1; % smoothness level

for k = 2:T
    f = randn(N,1);
    f = f / norm(f) * epsilon; 
    fdc = randn * 0 * ones(N,1); % DC component
    Temp(:,k) = Temp(:,k-1) + LHalfInv*f + fdc;
end


%% sampling matrix
SampleNum = floor(N*0.4); % the number of sampled points at each time
SampleMatrix = zeros(N,T);
for i = 1:T
    SampleMatrix(randperm(N, SampleNum),i) = 1;
end 

noise = 0.1 * randn(size(Temp)); % measurement noise
save paramAWD Position A W D L Temp noise SampleMatrix


%% generate data for Expreiment_A_5b
Tempall = zeros(N,T,7);
Tempall(:,:,4) = Temp;

[V,Lambda]=eig(L);
Lambda(1,1)=0;
lambda = diag(Lambda);
lambdaHalfInv = 1./sqrt(lambda);
lambdaHalfInv(1) = 0;
LHalfInv = V*diag(lambdaHalfInv)*V';

epsilon_set = [0.1 0.2 0.5 1 2 5 10]; % smoothness level
for i_epsilon = [1 2 3 5 6 7]
    epsilon = epsilon_set(i_epsilon);
    Tempall(:,1,i_epsilon) = Tempall(:,1,4);

    for k = 2:T
        f = randn(N,1);
        f = f / norm(f) * epsilon;
        fdc = randn * 0 * ones(N,1); % DC component
        Tempall(:,k,i_epsilon) = Tempall(:,k-1,i_epsilon) + LHalfInv*f + fdc;
    end    
end

save paramAWDall Position A W D L Tempall noise
    