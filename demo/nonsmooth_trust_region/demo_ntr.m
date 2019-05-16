clear all;
close all;
clc;

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset','*_circle_original.png','*_circle_noisy.png');

% Define lower level problem
lower_level_problem.solve = @(lambda,dataset) solve_ntr_lower_level(lambda,dataset.get_corrupt(1));

% Define upper level problem
upper_level_problem.eval = @(u,lambda,dataset) eval_ntr_upper_level(u,lambda,dataset);
upper_level_problem.gradient = @(u,lambda,dataset,params) solve_ntr_gradient(u,lambda,dataset,params);
upper_level_problem.dataset = dataset;

%% Solving the bilevel problem
bilevel_param.verbose = 2;
bilevel_param.maxit = 100;
bilevel_param.tol = 1e-3;
bilevel_param.algo = 'NONSMOOTH_TRUST_REGION';
bilevel_param.radius = 2.0;
bilevel_param.minradius = 0.1;
bilevel_param.gamma1 = 0.5;
bilevel_param.gamma2 = 2.0;
bilevel_param.eta1 = 0.10;
bilevel_param.eta2 = 0.80;
bilevel_param.use_bfgs = true;
%bilevel_param.use_sr1 = true;
lambda = 1;
[optimal_parameter,info] = solve_bilevel(lambda,lower_level_problem,upper_level_problem,bilevel_param);

%% Plotting
imagesc_gray(dataset.get_target(1),1,'Original','131');
imagesc_gray(dataset.get_corrupt(1),1,'Gaussian Noise','132');
imagesc_gray(info.u_history(:,:,end),1,'ROF Optimal Image Denoising','133');

%% Saving experiment
%save_experiment(info,'/Users/dvillacis/Desktop/ntr_experiment_sr1');
