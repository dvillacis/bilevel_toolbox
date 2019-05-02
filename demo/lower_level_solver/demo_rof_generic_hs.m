clear all;
clc;

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset','*_circle_original.png','*_circle_noisy.png');

%% Load image
noisy = dataset.get_corrupt(1);
[M,N] = size(noisy);

%% Solver Parameters
param.verbose = 2;
param.maxiter = 2000;
param.check = 200;

id_op = IdentityOperator([M,N]);
grad_op = FinDiffOperator([M,N],'fn');
q = zeros(M,N,2);

alpha = 1;
lambda = 4.5;
gamma = 1000;

%% Solving using Hintermuller-Stadler
[denoised,gap] = solve_generic_l1_l2_hs({lambda},{alpha},{id_op},{grad_op},{noisy},{q},{gamma},0*noisy,param);

%% Plotting
imagesc_gray(dataset.get_target(1),1,'Original','131');
imagesc_gray(noisy,1,'Gaussian Noise','132');
imagesc_gray(denoised,1,'ROF Image Denoising','133');