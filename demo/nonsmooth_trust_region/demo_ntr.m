clear all;
clc;

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset_single_gaussian','*_circle_original.png','*_circle_noisy.png');

% Define lower level problem
lower_level_problem.solve = @(lambda,dataset) solve_ntr_lower_level(lambda,dataset.get_corrupt(1));

% Define upper level problem
upper_level_problem.eval = @(u,lambda,dataset) eval_ntr_upper_level(u,lambda,dataset);
upper_level_problem.gradient = @(u,lambda,dataset,params) solve_ntr_gradient(u,lambda,dataset,params);
upper_level_problem.dataset = dataset;

%% Solving the bilevel problem
bilevel_param.verbose = 2;
bilevel_param.maxit = 400;
bilevel_param.tol = 1e-3;
bilevel_param.algo = 'NONSMOOTH_TRUST_REGION';
bilevel_param.radius = 4.0;
bilevel_param.minradius = 1.0;
bilevel_param.gamma1 = 0.5;
bilevel_param.gamma2 = 1.5;
bilevel_param.eta1 = 0.10;
bilevel_param.eta2 = 0.80;
bilevel_param.use_bfgs = true;
lambda = 15;
[optimal_parameter,info] = solve_bilevel(lambda,lower_level_problem,upper_level_problem,bilevel_param);

%% Saving experiment
%save_experiment(info,'/Users/dvillacis/Desktop/ntr_experiment');