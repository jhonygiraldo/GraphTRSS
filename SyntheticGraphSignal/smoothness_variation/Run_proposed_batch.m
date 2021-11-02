% ========================================================================
% Time-varying Graph Signal Reconstruction,
%
% Copyright(c) 2017 Kai Qiu and Yuantao Gu
% All Rights Reserved.
% ----------------------------------------------------------------------
% Batch Reconstruction of Time-Varying Graph Signal (BR-TVGS), objective function:
% f(x)= 0.5* ||J*x-y||^2 + 0.5*alpha* (Tx)'*L*(Tx) + 0.5*beta* x'*L*x + 0.5*gamma*||T*x||^2,
% where beta and gamma equal 0.
% 
% Version 1.0
% Written by Kai Qiu (q1987k@163.com)
%----------------------------------------------------------------------


clear; clc; close all;
addpath('../../../solvers');
addpath('../../../utilities');
load ../paramAWD
[N,T,M] = size(Tempall);

test_time = 50; % the sampling matrix is different at each time
epsilon_set = [0.1 0.2 0.5 1 2 5 10]; % smoothness level
for i_epsilon = 1 : length(epsilon_set)
    i_epsilon
    Temp = Tempall(:,:,i_epsilon);
    i_test = 1;
    SampleNum = floor(N*0.4); % the number of sampled points at each time
    SampleMatrix = zeros(N,T); % sampling matrix
    for i = 1:T
        SampleMatrix(randperm(N, SampleNum),i) = 1;
    end

    param.J = SampleMatrix(:,1:T); % sampling matrix
    param.y = param.J .* (Temp+noise); % sampled data
    param.L = L; % Graph Laplacian
    param.T = TV_Temp(); % temporal difference operator
    param.beta = 0; 
    param.gamma = 0; 
    param.niter = 10000;

    alpha_set = [1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 1e2, 2e2, 5e2];
    errorSeekRMSE = zeros(1, length(alpha_set));

    for i_alpha = 1 : length(alpha_set)      
        param.alpha = alpha_set(i_alpha);            
        param.x0 = 0 * param.y;
        x_recon = solver_BR_TVGS(param);
        errorSeekRMSE(i_alpha) = norm(Temp(:) - x_recon(:))/sqrt(N*T); %RMSE                      
    end
    [Error_RMSE(i_test,i_epsilon), Indx] = min(errorSeekRMSE);  
    
    
    for i_test = 2:test_time
        i_test
        SampleMatrix = zeros(N,T);
        for i = 1:T
            SampleMatrix(randperm(N, SampleNum),i) = 1;
        end

        param.J = SampleMatrix(:,1:T);
        param.y = param.J .* (Temp+noise);     
        param.alpha = alpha_set(Indx);
        param.x0 = 0 * param.y;
        x_recon = solver_BR_TVGS(param);
        Error_RMSE(i_test,i_epsilon) = norm(Temp(:) - x_recon(:))/sqrt(N*T); %RMSE       
    end
end

save Error_batch_RMSE Error_RMSE epsilon_set
