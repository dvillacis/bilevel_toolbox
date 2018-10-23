clear all;
close all;
clc;

init_bilevel_toolbox();

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset_single_gaussian','*_circle_original.png','*_circle_noisy.png');

%% Load input image
original = dataset.get_target(1);
noisy = dataset.get_corrupt(1);

% Define lower level problem
lower_level_problem.solve = @(alpha) solve_lower_level(alpha,noisy);

% Define upper level problem
upper_level_problem.eval = @(u,alpha) 0.5*norm(u-original(:)).^2;
upper_level_problem.adjoint = @(u,alpha,radius) solve_adjoint(u,alpha,radius);
upper_level_problem.slack = @(u,alpha) solve_slack(u,alpha,noisy);

%% Solving the bilevel problem
bilevel_param.verbose = 2;
bilevel_param.maxit = 100;
bilevel_param.tol = 1e-3;
bilevel_param.algo = 'NONSMOOTH_TRUST_REGION';
bilevel_param.radius = 0.1;
bilevel_param.gamma1 = 0.5;
bilevel_param.gamma2 = 1.5;
bilevel_param.eta1 = 0.25;
bilevel_param.eta2 = 0.75;
alpha = 0.5;
[sol,info] = solve_bilevel(alpha,lower_level_problem,upper_level_problem,bilevel_param);

function y = solve_lower_level(alpha,noisy)
  param_lower_level.maxit = 2000;
  param_lower_level.alpha = alpha;
  y = solve_rof_cp_single_gaussian(noisy,param_lower_level);
end

function q = solve_slack(u,alpha,noisy)
    [m,n] = size(noisy);
    A = gradient_matrix(m,n);
    q = A'\(noisy(:)-u(:));
    q = reshape(q,m*n,2);
end

function grad = solve_adjoint(u,alpha,radius)

end
