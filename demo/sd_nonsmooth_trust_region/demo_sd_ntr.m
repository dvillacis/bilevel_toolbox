clear all;
close all;
clc;

init_bilevel_toolbox();

%% Load dataset
%dataset = DatasetInFolder('data/circle_dataset','*_circle_original.png','*_circle_noisy.png');
dataset = DatasetInFolder('data/smiley','*_smiley_original.png','*_smiley_noisy.png');

%% Load input image
% original = dataset.get_target(1);
% noisy = dataset.get_corrupt(1);
[M,N]=size(dataset.get_target(1));

% Define lower level problem
lower_level_problem.solve = @(lambda,dataset) solve_sd_ntr_lower_level(lambda,dataset.get_corrupt(1));

% Define upper level problem
upper_level_problem.eval = @(u,lambda,dataset) eval_sd_ntr_upper_level(u,lambda,dataset);
upper_level_problem.gradient = @(u,lambda,dataset,params) solve_sd_ntr_gradient(u,lambda,dataset,params);
upper_level_problem.dataset = dataset;

%% Solving the bilevel problem
bilevel_param.verbose = 2;
bilevel_param.maxit = 400;
bilevel_param.tol = 1e-2;
bilevel_param.algo = 'NONSMOOTH_TRUST_REGION';
bilevel_param.radius = 1.0;
bilevel_param.minradius = 0.1;
bilevel_param.gamma1 = 0.5;
bilevel_param.gamma2 = 2.0;
bilevel_param.eta1 = 0.10;
bilevel_param.eta2 = 0.80;

%lambda = 80.0*triu(ones(M,N))+0.9*tril(ones(M,N)); % Initial guess
lambda = 2*ones(M,N);
[sol,info] = solve_bilevel(lambda,lower_level_problem,upper_level_problem,bilevel_param);

%% Plot Surface
figure(1)
[a,b] = meshgrid(1:M,1:N);
surf(a,b,sol);
matlab2tikz('/Users/dvillacis/OneDrive - Escuela Politécnica Nacional/DOCTORADO/ARTICLES/Bilevel_SD_ROF/plots/sd_parameter.tex');

%% Plot Images
imagesc_gray(dataset.get_target(1),2,'Original','221');
imagesc_gray(dataset.get_corrupt(1),2,'Gaussian Noise','222');
imagesc_gray(info.sol_history(:,:,end),2,'Learned Parameter','223');
imagesc_gray(info.u_history(:,:,end),2,'ROF Optimal Image Denoising','224');
imagesc_write(info.u_history(:,:,end),'/Users/dvillacis/OneDrive - Escuela Politécnica Nacional/DOCTORADO/ARTICLES/Bilevel_SD_ROF/plots/sd_parameter.tex'
