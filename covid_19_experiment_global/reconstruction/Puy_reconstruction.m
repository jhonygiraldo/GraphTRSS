clear all, close all, clc;
%%
load('../graph_construction/full_graph.mat');
load('../covid_19_new_cases.mat');
x_matrix = Data;
%%
load('../sampling/puy_pattern.mat');
signals_t = size(x_matrix,2);
%% Bandwidth
G = gsp_estimate_lmax(G);
% Fourier basis
G = gsp_compute_fourier_basis(G);
x_hat_square = (G.U'*x_matrix(:,end)).^2;
integral_x_hat_square = cumsum(x_hat_square);
band_limit = find(integral_x_hat_square <= 0.9*integral_x_hat_square(end));
band_limit = band_limit(end);
Uk = G.U(:,1:band_limit);
weights = sum(Uk.^2, 2)/sum(Uk(:).^2);
P = sparse(1:G.N, 1:G.N, 1./sqrt(weights), G.N, G.N);
%% lambda and epsilon optimization
reg_set = [1e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 1e2, 2e2, 5e2];
error_sobolev_batch_puy_opt_lambda = cell(1,length(reg_set));
%%
repetitions = 1;
for ii=1:repetitions
    ii
    for i=1:length(m)
        SampleMatrix = puy_patterns_opt_lambda{ii,i};
        x_recon = zeros(G.N,signals_t);
        for jj=1:length(reg_set)
            for j=1:signals_t
                ind_obs = find(SampleMatrix(:,j) == 1);
                nb_meas = length(ind_obs);
                M = sparse(1:nb_meas, ind_obs, 1, nb_meas, G.N);
                ynoise_init = M*x_matrix(:,j);
                % Prepare matrices for reconstruction
                MP = M*P; MtM = MP'*MP;
                ynoise = P*P*M'*ynoise_init;
                x_recon(:,j) = (MtM + reg_set(jj)*G.L)\ynoise;
            end
            indx_non_sampled = find(SampleMatrix(:) == 0);
            x_vector_original = x_matrix(indx_non_sampled);
            x_vector_reconstructed = x_recon(indx_non_sampled);
            error_sobolev_batch_puy_opt_lambda{1,jj}(ii,i) = mean(abs(x_vector_original-x_vector_reconstructed).^2);
        end
    end
end
%% Best parameters
mean_errors_sobolev_batch = zeros(size(error_sobolev_batch_puy_opt_lambda));
for i=1:size(error_sobolev_batch_puy_opt_lambda,1)
    for j=1:size(error_sobolev_batch_puy_opt_lambda,2)
        mean_errors_sobolev_batch(i,j) = mean(mean(error_sobolev_batch_puy_opt_lambda{i,j}));
    end
end
[~, best_reg] = find(mean_errors_sobolev_batch == min(min(mean_errors_sobolev_batch)));
%% Experiment best parameters
repetitions = 100;
RMSE = zeros(repetitions,length(m));
MAE = zeros(repetitions,length(m));
MAPE = zeros(repetitions,length(m));
reg_param = reg_set(best_reg);
for ii=1:repetitions
    ii
    for i=1:length(m)
        SampleMatrix = puy_patterns_experiment{ii,i};
        x_recon = zeros(G.N,signals_t);
        for j=1:signals_t
            ind_obs = find(SampleMatrix(:,j) == 1);
            nb_meas = length(ind_obs);
            M = sparse(1:nb_meas, ind_obs, 1, nb_meas, G.N);
            ynoise_init = M*x_matrix(:,j);
            % Prepare matrices for reconstruction
            MP = M*P; MtM = MP'*MP;
            ynoise = P*P*M'*ynoise_init;
            x_recon(:,j) = (MtM + reg_set(jj)*G.L)\ynoise;
        end
        %%
        indx_non_sampled = find(SampleMatrix(:) == 0);
        x_vector_original = x_matrix(indx_non_sampled);
        x_vector_reconstructed = x_recon(indx_non_sampled);
        RMSE(ii,i) = sqrt(mean((x_vector_original-x_vector_reconstructed).^2));
        MAE(ii,i) = mean(abs(x_vector_original-x_vector_reconstructed));
        index_zero = find(x_vector_original == 0);
        x_vector_original(index_zero) = [];
        x_vector_reconstructed(index_zero) = [];
        MAPE(ii,i) = mean(abs((x_vector_original-x_vector_reconstructed)./x_vector_original));
    end
end
%%
results_path = '../results/';
mkdir(results_path);
save([results_path 'error_puy.mat'],'error_sobolev_batch_puy_opt_lambda',...
    'RMSE','MAE','MAPE','m','reg_set','best_reg');