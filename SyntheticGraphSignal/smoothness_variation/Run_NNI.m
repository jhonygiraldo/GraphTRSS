% comparison method: Natural neighbor interpolation

clear; clc; close all;
addpath('../../../solvers');
load ../paramAWDall
[N,T,M] = size(Tempall);

test_time = 50; % the sampling matrix is different at each time
epsilon_set = [0.1 0.2 0.5 1 2 5 10]; % smoothness level
for i_epsilon = 1 : length(epsilon_set)
    i_epsilon
    Temp = Tempall(:,:,i_epsilon);   
    for i_test = 1 : test_time
        i_test       
        SampleNum = floor(N*0.4); % the number of sampled points at each time
        SampleMatrix = zeros(N,T); % sampling matrix
        for i = 1:T
            SampleMatrix(randperm(N, SampleNum),i) = 1;
        end

        TempS = SampleMatrix(:,1:T) .* (Temp+noise); % sampled data

        %% interpolated by function scatteredInterpolant
        x_recon = solver_NNI(SampleMatrix, TempS, Position);

        Error_RMSE(i_test,i_epsilon) = norm(Temp(:) - x_recon(:))/sqrt(N*T); %RMSE
    end
end

save Error_NNI_RMSE Error_RMSE epsilon_set
