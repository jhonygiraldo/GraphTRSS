clear all, close all, clc;
%%
load ../paramAWD
x_matrix = Temp;
%%
m = [0.1 : 0.1 : 0.9];  %Sampling density
N = size(x_matrix,1);
T = size(x_matrix,2);
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
x_hat_square = (G.U'*x_matrix(:,end)).^2;
integral_x_hat_square = cumsum(x_hat_square);
band_limit = find(integral_x_hat_square <= 0.9*integral_x_hat_square(end));
band_limit = band_limit(end);
Uk = G.U(:,1:band_limit);
weights = sum(Uk.^2, 2)/sum(Uk(:).^2);
%% Sampling patterns for experiments
repetitions = 50;
puy_patterns_experiment = cell(repetitions,length(m));
for ii=1:repetitions
    for i=1:length(m)
        random_pattern_temp = zeros(N,T);
        matrix_densities_truncated = weights;
        probability_distribution = matrix_densities_truncated/(sum(matrix_densities_truncated));
        num_samples = round(N*m(i));
        for jj=1:T
            indx_sampled = datasample(1:N,num_samples,'Replace',false,'Weights',probability_distribution);
            random_pattern_temp(indx_sampled,jj) = 1;
        end
        puy_patterns_experiment{ii,i} = random_pattern_temp;
    end
end
%%
save('puy_pattern.mat','puy_patterns_experiment','m');