% ========================================================================
% Time-varying Graph Signal Reconstruction,
%
% Copyright(c) 2017 Kai Qiu and Yuantao Gu
% All Rights Reserved.
% ----------------------------------------------------------------------
% Online Reconstruction of Time-Varying Graph Signal (OR-TVGS), objective function:
% f(x) = 0.5 * ||J*x - y||^2 + 0.5 * alpha * (x-xref)'*L*(x-xref).
% 
% Version 1.0
% Written by Kai Qiu (q1987k@163.com)
%----------------------------------------------------------------------


clear; clc; close all;
addpath('../../../solvers');
load ../paramAWD
[N,T] = size(Temp);

noise_set = [0.01 0.1 0.2 0.4 0.6 0.8 1]; 
for i_noise = 1 : length(noise_set)
    i_noise
    noise = noise_set(i_noise) * randn(size(Temp));   
    y_all = SampleMatrix .* (Temp+noise); % sampled data
    param1.L = L; % Graph Laplacian
    param1.niter = 10000;

    alpha_set = [1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50];
    errorSeekRMSE = zeros(1, length(alpha_set));
    for i_alpha = 1 : length(alpha_set)   
        i_alpha
        param1.alpha = alpha_set(i_alpha);
        param1.ref = zeros(N,1);
        
        x_recon = zeros(size(Temp));
        for i=1:T
            param1.J = SampleMatrix(:,i);
            param1.y = y_all(:,i);
            
            param1.x0 = param1.y; % the reference signal is the sampled signal for i=1.
            if i >1
                % the reference signal is the reconstructed signal of the previous time for i>1.
                param1.x0 = x_recon(:,i-1); 
            end
            x_recon(:,i) = solver_OR_TVGS(param1);
            param1.ref = x_recon(:,i);
        end
        errorSeekRMSE(i_alpha) = norm(Temp(:) - x_recon(:))/sqrt(N*T); %RMSE
    end
    Error_RMSE(i_noise) = min(errorSeekRMSE);  
end

save Error_online_RMSE Error_RMSE noise_set
