function opts = get_opts_FPCA(xs,maxr,m,n,sr,fr)

% Shiqian Ma
% Apr. 24, 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Here is a simple way to determine whether this is an easy problem
%%% By "easy", we mean that the rank of the matrix is very small compared 
%%% to the dimension and number of samples, so the problem can be solved by
%%% SDP solvers. In this case, all parameters can be set very loose to
%%% achieve a very fast speed. If it does not work, which means this is
%%% actually a very "hard" problem, where compared to the dimension and
%%% number of samples, the rank of the matrix is not relatively large 
%%% (usually, such problem cannot be solved by SDP solvers), then we
%%% should set some strict parameters to make sure that we can solve the
%%% problem correctly. 
%%%
%%% In summary, the parameters for "hard = 0" guarantee that the solver is
%%% fast; the parameters for "hard = 1" guarantee that it can get great
%%% recoverability. The users can switch between "hard = 0" and "hard = 1"
%%% and change the parameters if necessary. 
%%%
if sr <= 0.5*fr || fr >= 0.38
    hard = 1;
%     fprintf('This is a "hard" problem! \n\n');
else
    hard = 0;
%     fprintf('This is an "easy" problem! \n\n');
end
hard = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if hard && max(m,n) < 1000
    opts.mu = 1e-8; % final mu
    opts.xtol = 1e-8; % tolerance for subproblems in continuation
    opts.maxinneriter = 500; % maximum iteration number for subproblems in continuation
    opts.tau = 1; % stepsize
else
    opts.mu = 1e-4;
    opts.xtol = 1e-4;
    opts.maxinneriter = 10;
    opts.tau = 2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
opts.xs = xs;  % true solution xs is given
opts.fastsvd_ratio_leading = 1e-2; % ratio for computing the hard thresholding, i.e., the rank for next iteration
opts.mxitr = 1e+5; % maximum iteration number for the outer loop
opts.eta = 1/4; % ratio to decrease mu in continuation
opts.maxr = maxr; % maximum rank r is given 
opts.print = 0;
