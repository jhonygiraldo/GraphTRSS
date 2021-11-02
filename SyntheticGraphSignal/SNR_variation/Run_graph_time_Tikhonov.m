% comparison method: Graph-time Tikhonov
% Objective function:
% f(x)= 0.5* ||J*x-y||^2 + 0.5*alpha* (Tx)'*L*(Tx) + 0.5*beta* x'*L*x + 0.5*gamma*||T*x||^2,
% where alpha equals 0.

clear; clc; close all;
addpath('../../../solvers');
addpath('../../../utilities');
load ../paramAWD
[N,T] = size(Temp);

noise_set = [0.01 0.1 0.2 0.4 0.6 0.8 1];
for i_noise = 1 : length(noise_set)
    i_noise
    noise = noise_set(i_noise) * randn(size(Temp));

    param.J = SampleMatrix(:,1:T); % sampling matrix
    param.y = param.J .* (Temp+noise); % sampled data
    param.L = L; % Graph Laplacian
    param.T = TV_Temp(); % temporal difference operator
    param.alpha = 0;
    param.niter = 10000;

    beta_set = [1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5];
    gamma_set = [1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5];
    errorSeekRMSE = zeros(length(beta_set), length(gamma_set));

    for i_beta = 1 : length(beta_set)
        param.beta = beta_set(i_beta);
        for i_gamma = 1 : length(gamma_set)
            param.gamma = gamma_set(i_gamma);           
            param.x0 = 0 * param.y;
            x_recon = solver_BR_TVGS(param);
            errorSeekRMSE(i_beta,i_gamma) = norm(Temp(:) - x_recon(:))/sqrt(N*T); %RMSE
        end
    end
    Error_RMSE(i_noise) = min(errorSeekRMSE(:));  
end

save Error_station_RMSE Error_RMSE noise_set
