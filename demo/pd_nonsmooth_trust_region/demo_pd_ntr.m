clear all;
close all;
clc;

init_bilevel_toolbox();

%% Load dataset
%dataset = DatasetInFolder('data/circle_dataset_single_gaussian','*_circle_original.png','*_circle_noisy.png');
dataset = DatasetInFolder('data/smiley','*_smiley_original.png','*_smiley_noisy.png');

%% Load input image
% original = dataset.get_target(1);
% noisy = dataset.get_corrupt(1);
[M,N]=size(dataset.get_target(1));

% Define lower level problem
lower_level_problem.solve = @(lambda,dataset) solve_pd_ntr_lower_level(lambda,dataset.get_corrupt(1));

% Define upper level problem
upper_level_problem.eval = @(u,lambda,dataset) eval_pd_ntr_upper_level(u,lambda,dataset);
upper_level_problem.gradient = @(u,lambda,dataset,params) solve_pd_ntr_gradient(u,lambda,dataset,params);
upper_level_problem.dataset = dataset;

%% Solving the bilevel problem
bilevel_param.verbose = 2;
bilevel_param.maxit = 300;
bilevel_param.tol = 1e-4;
bilevel_param.algo = 'NONSMOOTH_TRUST_REGION';
bilevel_param.radius = 500.0;
bilevel_param.minradius = 10.0;
bilevel_param.gamma1 = 0.5;
bilevel_param.gamma2 = 1.5;
bilevel_param.eta1 = 0.10;
bilevel_param.eta2 = 0.60;
bilevel_param.use_bfgs = true;
%lambda = 80.0*triu(ones(M,N))+0.9*tril(ones(M,N)); % Initial guess
lambda = 100*ones(32,32);
[sol,info] = solve_bilevel(lambda,lower_level_problem,upper_level_problem,bilevel_param);
po = PatchOperator(size(sol),[M,N]);

%% Plot
figure(1)
[a,b] = meshgrid(1:M,1:N);
surf(a,b,po.val(sol));

figure(2)
imagesc_gray(info.u_history(:,:,end));

figure(3)
imagesc_gray(po.val(info.sol_history(:,:,end)));

%% Save animations
%create_animation_image(info.u_history,'/Users/dvillacis/Desktop/smiley_pd_evolution2.gif');
%create_animation_parameter(info.sol_history,'/Users/dvillacis/Desktop/smiley_pd_param_evolution2.gif',po);


