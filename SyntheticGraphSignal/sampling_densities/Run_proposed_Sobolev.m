clear; clc; close all;
addpath('../../../solvers');
addpath('../../../utilities');
load ../paramAWD
[N,T] = size(Temp);

%% alpha parameter
alpha_set = [1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 1e2, 2e2, 5e2];
epsilon_set = [1e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 1e2, 2e2, 5e2];
sobolev_term_matrices = cell(1,length(epsilon_set));
for i=1:length(epsilon_set)
    laplacian = L+epsilon_set(i)*eye(N);
    %% Symmetrization
    sobolev_term = 0.5*(laplacian+laplacian.');
    sobolev_term_matrices{1,i} = sparse(sobolev_term);
end

test_time = 50; % the sampling matrix is different at each time
rate_set = 0.1 : 0.1 : 0.9; % sampling rate
best_epsilon = zeros(1,length(rate_set));
best_alpha = zeros(1,length(rate_set));
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
    %param.L = L; % Graph Laplacian
    param.T = TV_Temp(); % temporal difference operator
    param.beta = 0; 
    param.gamma = 0; 
    param.niter = 10000;
    error_sobolev_batch_random = zeros(length(epsilon_set),length(alpha_set));
    %%
    for i_epsilon=1 : length(epsilon_set)
        param.L = sobolev_term_matrices{1,i_epsilon};
        for i_alpha = 1 : length(alpha_set)      
            param.alpha = alpha_set(i_alpha);            
            param.x0 = 0 * param.y;
            x_recon = solver_BR_TVGS(param);
            error_sobolev_batch_random(i_epsilon,i_alpha) = norm(Temp(:) - x_recon(:))/sqrt(N*T); %RMSE                      
        end
    end
    [Indx_Epsilon, Indx_alpha] = find(error_sobolev_batch_random == min(min(error_sobolev_batch_random)));
    best_epsilon(i_rate) = epsilon_set(Indx_Epsilon);
    best_alpha(i_rate) = alpha_set(Indx_alpha);
    Error_RMSE(i_test,i_rate) = min(min(error_sobolev_batch_random));
    
    for i_test = 2:test_time
        i_test
        SampleMatrix = zeros(N,T);
        for i = 1:T
            SampleMatrix(randperm(N, SampleNum),i) = 1;
        end

        param.J = SampleMatrix(:,1:T);
        param.L = sobolev_term_matrices{1,Indx_Epsilon};
        param.y = param.J .* (Temp+noise);
        param.alpha = alpha_set(Indx_alpha);
        param.x0 = 0 * param.y;
        x_recon = solver_BR_TVGS(param);
        Error_RMSE(i_test,i_rate) = norm(Temp(:) - x_recon(:))/sqrt(N*T); %RMSE       
    end
end

save Error_Sobolev_RMSE Error_RMSE rate_set best_epsilon best_alpha