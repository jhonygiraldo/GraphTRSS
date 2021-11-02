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
[N,T] = size(Temp);

noise_set = [0.01 0.1 0.2 0.4 0.6 0.8 1]; 
for i_noise = 1 : length(noise_set)
    i_noise
    noise = noise_set(i_noise) * randn(size(Temp));
    
    param.J = SampleMatrix(:,1:T);
    param.y = param.J .* (Temp+noise);
    param.L = L;
    param.T = TV_Temp();
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
    Error_RMSE(i_noise) = min(errorSeekRMSE);  
end

save Error_batch_RMSE Error_RMSE noise_set
