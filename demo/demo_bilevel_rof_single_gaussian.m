clear all;
close all;
clc;

verbose = 2;

init_bilevel_toolbox();

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset_single_gaussian','*_circle_original.png','*_circle_noisy.png');

%% Load input image
original = dataset.get_target(1);
noisy = dataset.get_corrupt(1);

%% Setup for the lower level problem
lower_f2.grad = @(u) u - noisy;
lower_f2.eval = @(u) 0.5 * norm(u-noisy)^2;
lower_f2.beta = 1;

lower_param_tv.verbose = verbose - 1;
lower_alpha = 0.1;
lower_f1.prox = @(u,T) prox_tv(u,alpha*T,lower_param_tv);
lower_f1.eval = @(u) alpha*norm_tv(u);

% Lower level parameters
lower_param_solver.verbose = verbose;
lower_param_solver.maxit = 300;
lower_param_solver.tol = 1e-7;
lower_param_solver.gamma = 0.1;

% Setting up the lower level problem data structure
lower_level_problem.f2 = lower_f2;
lower_level_problem.f1 = lower_f1;
lower_level_problem.alpha = lower_alpha;
lower_level_problem.param_tv = lower_param_tv;
lower_level_problem.param = lower_param_solver;

%% Setup for the upper level problem
upper_f1.eval = @(u) 0.5 * norm(u-original)^2;
upper_f1.grad = @(u) u - original;
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
bilevel_param.tol = 1e-7;
bilevel_param.algo = 'NONSMOOTH_TRUST_REGION';
[sol,info] = solve_bilevel(lower_alpha,lower_level_problem,upper_level_problem,bilevel_param);
