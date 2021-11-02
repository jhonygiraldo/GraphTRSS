function X_recon = solver_batch_sobolev_min(param)
% This function solves the reconstruction of time-varying graph signals
% using Sobolev norm minimization using gradient prohection algorithm.
X_previous = param.X_0;
for i=1:param.number_iter
    gradient_i = param.L*X_previous*param.D_h*param.D_h';
    step_tau = param.tau;
    X_current = X_previous - step_tau*gradient_i; % Gradient descent step
    X_current = param.Y+X_current-(param.J .* X_current); % Projection of answer
    X_current = max(0,X_current); % Projection of answer
    while trace((X_current*param.D_h)'*param.L*X_current*param.D_h) > ...
            trace((X_previous*param.D_h)'*param.L*X_previous*param.D_h) + ...
            param.alpha*trace(gradient_i'*(X_current-X_previous)) & i>1
        %dot(gradient_i,X_current-X_previous)
        step_tau = param.beta*step_tau;
        X_current = X_previous-step_tau*gradient_i; % Gradient descent step
        X_current = param.Y+X_current-(param.J .* X_current); % Projection of answer
        X_current = max(0,X_current); % Projection of answer
        if norm(X_current-X_previous) < 0.05
            break;
        end
    end
    %disp(['Difference of norms ',num2str(norm(X_current-X_previous))]);
    if norm(X_current-X_previous) < 0.05
        break;
    end
    X_previous = X_current;
    %disp(['Objective function ',num2str(trace((X_current*param.D_h)'*param.L*X_current*param.D_h))]);
end
disp(['Iterations: ',num2str(i)]);
X_recon = X_current;