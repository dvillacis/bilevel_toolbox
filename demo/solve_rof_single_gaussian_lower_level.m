function [sol,info] = solve_rof_single_gaussian_lower_level(noisy,alpha)
%SOLVE_ROF_SINGLE_GAUSSIAN_LOWER_LEVEL Summary of this function goes here
%   Detailed explanation goes here

    lower_f2.grad = @(u) u - noisy;
    lower_f2.eval = @(u) 0.5 * norm(u-noisy)^2;
    lower_f2.beta = 1;

    lower_param_tv.verbose = 0;
    lower_f1.alpha = alpha;
    lower_f1.prox = @(u,T) prox_tv(u,alpha*T,lower_param_tv);
    lower_f1.eval = @(u) alpha*norm_tv(u);

    % Lower level parameters
    lower_param_solver.verbose = 0;
    lower_param_solver.maxit = 300;
    lower_param_solver.tol = 1e-8;
    lower_param_solver.gamma = 0.1;

    [sol,info] = solvep(zeros(size(noisy)),{lower_f1,lower_f2},lower_param_solver);
end

