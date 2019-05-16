clear all;
close all;
clc;

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset','*_circle_original.png','*_circle_noisy.png');

%% Load image
noisy = dataset.get_corrupt(4);
[M,N] = size(noisy);

%% Solver Parameters
param.verbose = 2;
param.maxiter = 2000;
param.check = 500;

id_op = IdentityOperator([M,N]);
grad_op = FinDiffOperator([M,N],'fn');
q = zeros(M,N,2);

alpha = ones(M,N);
lambda = [3 ; 1.2];
po = PatchOperator(size(lambda),[M,N]);
lambda_out = po.val(lambda);

%% Solving
[denoised,gap] = solve_generic_l1_l2({lambda_out},{alpha},{id_op},{grad_op},noisy,{q},0,0*noisy,param);

%% Plotting
imagesc_gray(dataset.get_target(4),1,'Original','141');
imagesc_gray(noisy,1,'Gaussian Noise','142');
imagesc_gray(lambda_out,1,'Lambda','143');
imagesc_gray(denoised,1,'PD-ROF Image Denoising','144');