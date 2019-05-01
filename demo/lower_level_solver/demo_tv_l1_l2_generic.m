clear all;
close all;
clc;

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset','*_circle_original.png','*_circle_noisy.png');

%% Load image
noisy = dataset.get_corrupt(3);
[M,N] = size(noisy);

%% Solver Parameters
param.verbose = 2;
param.maxiter = 5000;
param.check = 200;
param.tol = 1e-3;

id_op = IdentityOperator([M,N]);
grad_op = FinDiffOperator([M,N],'fn');

q = zeros(M,N,2);

alpha_1 = 1;
lambda_1 = 0.5;
lambda_2 = 0.8;

%% Solving
[denoised,gap] = solve_generic_l1_l2({lambda_1},{lambda_2,alpha_1},{id_op},{id_op,grad_op},noisy,{noisy,q},0,0*noisy,param);

%% Plotting
figure
imagesc_gray(denoised);
