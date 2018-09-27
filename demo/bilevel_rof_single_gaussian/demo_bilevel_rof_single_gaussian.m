clear all;
close all;
clc;

verbose = 2;

init_bilevel_toolbox();

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset_single_gaussian','*_circle_original.png','*_circle_noisy.png');

%% Load input image
% original = dataset.get_target(1);
% noisy = dataset.get_corrupt(1);

%% Setup for the lower level problem
lower_level_problem.solve = @(noisy,alpha) solve_rof_single_gaussian_lower_level(noisy,alpha);

%% Setup for the upper level problem
upper_f1.eval = @(u,original) 0.5 * norm(u-original)^2;
upper_f1.grad = @(u,original) u - original;
upper_f1.beta = 1;

% Upper level parameters
upper_param_solver.verbose = verbose;
upper_param_solver.maxit = 300;
upper_param_solver.tol = 1e-7;

% Setting the upper level problem data structure
upper_level_problem.f1 = upper_f1;
upper_level_problem.param = upper_param_solver;

%% Solving the bilevel problem
bilevel_param.verbose = verbose;
bilevel_param.maxit = 100;
bilevel_param.tol = 1e-3;
bilevel_param.algo = 'NONSMOOTH_TRUST_REGION';
bilevel_param.radius = 0.1;
bilevel_param.gamma1 = 0.5;
bilevel_param.gamma2 = 1.5;
bilevel_param.eta1 = 0.25;
bilevel_param.eta2 = 0.75;
alpha = 0.5;
[sol,info] = solve_bilevel(alpha,dataset,lower_level_problem,upper_level_problem,bilevel_param);
