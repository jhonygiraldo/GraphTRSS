function x = solver_OR_TVGS(param)
%
% Online Reconstruction of Time-Varying Graph Signal (OR-TVGS), using conjugate gradient method.
% 
% Given the sampling model y = J*x, minimizes the following objective function:
%
% f(x) = 0.5 * ||J*x - y||^2 + 0.5 * alpha * (x-xref)'*L*(x-xref)
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
%    ref   -- the reference signal
%
%    alpha -- regularization parameter
%
%    niter -- maximum number of iterations
%
%  Output:
%
%      x  -- the reconstructed graph signal
% 
%  Written by Kai Qiu (q1987k@163.com)
%  date 6/20/2017


% starting point
x = param.x0;
gradToll = 1e-6 ;
k = 0;

g0 = grad(x,param);
dx = -g0;

while(1)
    % Optimal stepsize decision
    tmp = grad(dx,param)+param.y+ param.alpha * param.L * param.ref;
    t=-(g0(:)'*dx(:))/(tmp(:)'*dx(:));

    % update x
	x = (x + t*dx);
    
    %conjugate gradient calculation
	g1 = grad(x,param);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;
	
	% stopping criteria (to be improved)
	if (k > param.niter) || (norm(dx(:)) < gradToll), break;end
end
return;


function g = grad(x,param)%***********************************************

% part 1
Grad1 = param.J.*x-param.y;

% part 2
Grad2 = 0;
if param.alpha
    Grad2 = param.L * (x-param.ref);   
end

% composite gradient
g = Grad1 + param.alpha * Grad2;

