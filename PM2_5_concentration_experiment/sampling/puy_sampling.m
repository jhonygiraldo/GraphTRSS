clear all, close all, clc;
%%
load('../graph_construction/full_graph.mat');
load('../PM2_5_concentration.mat');
x_matrix = myDataPM;
x_matrix = x_matrix(:,1:220);
%%
good_data = (x_matrix>0); % indicating the valid data
%%
m = [0.1:0.05:0.45];  %Sampling density
N = size(x_matrix,1);
T = size(x_matrix,2);
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
%% Sampling patterns for lambda optimization
repetitions = 1;
puy_patterns_opt_lambda = cell(repetitions,length(m));
for ii=1:repetitions
    for i=1:length(m)
        random_pattern_temp = zeros(N,T);
        matrix_densities_truncated = weights .* good_data(:,1);
        probability_distribution = matrix_densities_truncated/(sum(matrix_densities_truncated));
        num_samples = round(N*m(i));
        for jj=1:T
            indx_sampled = datasample(1:N,num_samples,'Replace',false,'Weights',probability_distribution);
            random_pattern_temp(indx_sampled,jj) = 1;
            if jj<T
                matrix_densities_truncated = weights .* good_data(:,jj+1);
                probability_distribution = matrix_densities_truncated/(sum(matrix_densities_truncated));
            end
        end
        puy_patterns_opt_lambda{ii,i} = random_pattern_temp;
    end
end
%% Sampling patterns for experiments
repetitions = 100;
puy_patterns_experiment = cell(repetitions,length(m));
for ii=1:repetitions
    for i=1:length(m)
        random_pattern_temp = zeros(N,T);
        matrix_densities_truncated = weights .* good_data(:,1);
        probability_distribution = matrix_densities_truncated/(sum(matrix_densities_truncated));
        num_samples = round(N*m(i));
        for jj=1:T
            indx_sampled = datasample(1:N,num_samples,'Replace',false,'Weights',probability_distribution);
            random_pattern_temp(indx_sampled,jj) = 1;
            if jj<T
                matrix_densities_truncated = weights .* good_data(:,jj+1);
                probability_distribution = matrix_densities_truncated/(sum(matrix_densities_truncated));
            end
        end
        puy_patterns_experiment{ii,i} = random_pattern_temp;
    end
end
%%
real_sampling_densities = zeros(size(puy_patterns_experiment));
for i=1:size(puy_patterns_experiment,1)
    for j=1:size(puy_patterns_experiment,2)
        real_sampling_densities(i,j) = sum(puy_patterns_experiment{i,j}(:))/(N*T);
    end
end
m_real = mean(real_sampling_densities);
%%
save('puy_pattern.mat','puy_patterns_experiment','puy_patterns_opt_lambda',...
    'm','m_real');