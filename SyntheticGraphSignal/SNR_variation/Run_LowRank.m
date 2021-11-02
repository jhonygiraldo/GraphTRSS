% comparison method: Low-rank matrix completion

clear; clc; close all;
addpath('../../../solvers');
addpath('../../../utilities');
load ../paramAWD
[N,T] = size(Temp);

noise_set = [0.01 0.1 0.2 0.4 0.6 0.8 1];
for i_noise = 1 : length(noise_set) %仿真不同的采样率
    i_noise
    noise = noise_set(i_noise) * randn(size(Temp));
    TempS = SampleMatrix(:,1:T) .* (Temp+noise);

    %%
    n = N*T;
    k = 16; % estimated rank (the results are not senstive to k)
    dr = k*(N+T-k);
    m = sum(SampleMatrix(:));

    loc = find(SampleMatrix(:)~=0); % locations of sampled points
    y = TempS(loc); % sampled values
    sr = m/n; % sampling rate
    maxr = floor(((N+T)-sqrt((N+T)^2-4*m))/2);
    opts = get_opts_FPCA(TempS,maxr,N,T,sr,dr/m); % get the parameters
    Out = solver_FPCA_MatComp(N,T,loc,y,opts);
    x_recon = Out.x;

    Error_RMSE(i_noise) = norm(Temp(:) - x_recon(:))/sqrt(N*T); %RMSE
end

save Error_LowRank_RMSE Error_RMSE noise_set
