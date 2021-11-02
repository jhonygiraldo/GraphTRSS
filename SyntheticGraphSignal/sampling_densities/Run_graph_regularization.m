% comparison method: Graph regularization
% Objective function:
% f(x)= 0.5* ||J*x-y||^2 + 0.5*alpha* (Tx)'*L*(Tx) + 0.5*beta* x'*L*x + 0.5*gamma*||T*x||^2,
% where alpha and gamma equal 0.

clear; clc; close all;
addpath('../../../solvers');
addpath('../../../utilities');
load ../paramAWD
[N,T] = size(Temp);



test_time = 50; % the sampling matrix is different at each time
rate_set = 0.1 : 0.1 : 0.9; % sampling rate
for i_rate = 1 : length(rate_set)
    i_rate
    i_test = 1;
    SampleNum = floor(N*rate_set(i_rate)); % the number of sampled points at each time
    SampleMatrix = zeros(N,T); % sampling matrix
    for i = 1:T
        SampleMatrix(randperm(N, SampleNum),i) = 1;
    end

    param.J = SampleMatrix(:,1:T); % sampling matrix
    param.y = param.J .* (Temp+noise); % sampled data
    param.L = L; % Graph Laplacian
    param.T = TV_Temp(); % temporal difference operator
    param.alpha = 0;
    param.gamma = 0;
    param.niter = 10000;

    beta_set = [1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5];
    errorSeekRMSE = zeros(1, length(beta_set));

    for i_beta = 1 : length(beta_set)                 
        param.beta = beta_set(i_beta);
        param.x0 = 0 * param.y;
        x_recon = solver_BR_TVGS(param);
        errorSeekRMSE(i_beta) = norm(Temp(:) - x_recon(:))/sqrt(N*T); %RMSE                      
    end
    [Error_RMSE(i_test,i_rate), Indx] = min(errorSeekRMSE);  
   
        
    for i_test = 2:test_time
        i_test
        SampleMatrix = zeros(N,T);
        for i = 1:T
            SampleMatrix(randperm(N, SampleNum),i) = 1;
        end

        param.J = SampleMatrix(:,1:T);
        param.y = param.J .* (Temp+noise);
        param.beta = beta_set(Indx);
        param.x0 = 0 * param.y;
        x_recon = solver_BR_TVGS(param);
        Error_RMSE(i_test,i_rate) = norm(Temp(:) - x_recon(:))/sqrt(N*T); %RMSE
    end
end

save Error_gsmooth_RMSE Error_RMSE rate_set
