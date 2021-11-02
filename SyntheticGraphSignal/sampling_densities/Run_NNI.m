% comparison method: Natural neighbor interpolation

clear; clc; close all;
addpath('../../../solvers');
load ../paramAWD
[N,T] = size(Temp);


test_time = 50; % the sampling matrix is different at each time
rate_set = 0.1 : 0.1 : 0.9; % sampling rate
for i_rate = 1 : length(rate_set)
    i_rate
    for i_test = 1 : test_time
        i_test       
        SampleNum = floor(N*rate_set(i_rate)); % the number of sampled points at each time
        SampleMatrix = zeros(N,T); % sampling matrix
        for i = 1:T
            SampleMatrix(randperm(N, SampleNum),i) = 1;
        end

        TempS = SampleMatrix(:,1:T) .* (Temp+noise); % sampled data

        %% interpolated by function scatteredInterpolant
        x_recon = solver_NNI(SampleMatrix, TempS, Position);

        Error_RMSE(i_test,i_rate) = norm(Temp(:) - x_recon(:))/sqrt(N*T); %RMSE
    end
end

save Error_NNI_RMSE Error_RMSE rate_set
