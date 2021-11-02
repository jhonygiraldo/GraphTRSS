clear; clc; close all;
load ../paramAWD
[N,T] = size(Temp);

load ../sampling_puy/puy_pattern

%% Create graph 
G.N = N;
G.W = W;
G.A = A;
G.coords = Position;
G.type = 'nearest neighbors';
G = gsp_graph_default_parameters(G);
%% Bandwidth
G = gsp_estimate_lmax(G);
% Fourier basis
G = gsp_compute_fourier_basis(G);
x_hat_square = (G.U'*Temp(:,end)).^2;
integral_x_hat_square = cumsum(x_hat_square);
band_limit = find(integral_x_hat_square <= 0.9*integral_x_hat_square(end));
band_limit = band_limit(end);
Uk = G.U(:,1:band_limit);
weights = sum(Uk.^2, 2)/sum(Uk(:).^2);
P = sparse(1:G.N, 1:G.N, 1./sqrt(weights), G.N, G.N);
%% lambda optimization
reg_set = [1e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 1e2, 2e2, 5e2];

SampleMatrix = puy_patterns_experiment{1,4};

noise_set = [0.01 0.1 0.2 0.4 0.6 0.8 1]; 
for i_noise = 1 : length(noise_set)
    i_noise
    noise = noise_set(i_noise) * randn(size(Temp));   
    y_all = SampleMatrix .* (Temp+noise); % sampled data
    
    errorSeekRMSE = zeros(1, length(reg_set));
    for i_reg_set = 1 : length(reg_set)   
        i_reg_set;
        
        x_recon = zeros(size(Temp));
        for j=1:T
            ind_obs = find(SampleMatrix(:,j) == 1);
            nb_meas = length(ind_obs);
            M = sparse(1:nb_meas, ind_obs, 1, nb_meas, G.N);
            ynoise_init = M*y_all(:,j);
            % Prepare matrices for reconstruction
            MP = M*P; MtM = MP'*MP;
            ynoise = P*P*M'*ynoise_init;
            x_recon(:,j) = (MtM + reg_set(i_reg_set)*G.L)\ynoise;
        end
        errorSeekRMSE(i_reg_set) = norm(Temp(:) - x_recon(:))/sqrt(N*T); %RMSE
    end
    Error_RMSE(i_noise) = min(errorSeekRMSE);  
end

save Error_puy_RMSE Error_RMSE noise_set