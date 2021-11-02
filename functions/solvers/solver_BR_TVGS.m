function [x, iterations, loss] = solver_BR_TVGS(param)
%
% Batch Reconstruction of Time-Varying Graph Signal (BR-TVGS), using conjugate gradient method.
% 
% Given the sampling model y = J*x, minimizes the following objective function:
%
% f(x)= 0.5* ||J*x-y||^2 + 0.5*alpha* (Tx)'*L*(Tx) + 0.5* beta* x'*L*x + 0.5*gamma* ||T*x||^2 
% 
%  Inputs:
%
%    x0  -- the initial signal for iterations
%   
%    J   -- the sampling operator
%
%    y   -- the sampled data
%
%    L   -- the Graph Laplacian
%
%    T   -- the temporal difference operator
%
%    alpha -- regularization parameter
%
%    beta  -- regularization parameter (default 0)
%
%    gamma -- regularization parameter (default 0)
%
%    niter -- maximum number of iterations
%
%  Output:
%
%      x  -- the reconstructed time-varying graph signal
% 
%  Written by Kai Qiu (q1987k@163.com)
%  date 6/20/2017


% starting point
x = param.x0;
gradToll = 1e-6 ;
k = 0;
g0 = grad(x,param);
dx = -g0;
loss = zeros(1,param.niter);

while(1)
    % Optimal stepsize decision
    tmp = grad(dx,param)+param.y;
    t=-(g0(:)'*dx(:))/(tmp(:)'*dx(:));

    % update x
	x = (x + t*dx);
    
    %conjugate gradient calculation
	g1 = grad(x,param);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;
    
    iterations = k;
	% stopping criteria (to be improved)
	if (k > param.niter) || (norm(dx(:)) < gradToll), 
        fprintf(' iteration number = %d\n',k);
        break;
    end
    loss(k) = norm(param.J.*x-param.y,'fro')^2 + (param.alpha/2)* trace((param.T * x)' * param.L * (param.T * x));
end
return;


function g = grad(x,param)%***********************************************

% part 1
Grad1 = param.J.*x-param.y;

% part 2
Grad2 = 0;
if param.alpha
    Grad2 =  param.T' * (param.L * (param.T * x));
end

% part 3
Grad3 = 0;
if param.beta
    Grad3 =  param.L * x;
end

% part 4
Grad4 = 0;
if param.gamma
    Grad4 =  param.T' * (param.T * x);
end

% composite gradient
g = Grad1 + param.alpha * Grad2 + param.beta * Grad3 + param.gamma * Grad4;