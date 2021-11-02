% comparison method: Natural neighbor interpolation

clear; clc; close all;
addpath('../../../solvers');
load ../paramAWD
[N,T] = size(Temp);

noise_set = [0.01 0.1 0.2 0.4 0.6 0.8 1];
for i_noise = 1 : length(noise_set)
    i_noise
    noise = noise_set(i_noise) * randn(size(Temp));
    TempS = SampleMatrix(:,1:T) .* (Temp+noise); % sampled data

    %% interpolated by function scatteredInterpolant
    x_recon = solver_NNI(SampleMatrix, TempS, Position);

    Error_RMSE(i_noise) = norm(Temp(:) - x_recon(:))/sqrt(N*T); %RMSE
end

save Error_NNI_RMSE Error_RMSE noise_set
