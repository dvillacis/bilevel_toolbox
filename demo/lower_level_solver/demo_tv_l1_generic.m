clear all;
close all;
clc;

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset','*_circle_original.png','*_circle_noisy.png');

%% Load image
noisy = dataset.get_corrupt(2);
[M,N] = size(noisy);

%% Solver Parameters
param.verbose = 2;
param.maxiter = 4000;
param.check = 500;

id_op = IdentityOperator([M,N]);
grad_op = FinDiffOperator([M,N],'fn');
q = zeros(M,N,2);

alpha = 1;
lambda = 0.8;

%% Solving
[denoised,gap] = solve_generic_l1_l2({},{lambda,alpha},{},{id_op,grad_op},{},{noisy,q},0,0*noisy,param);

%% Plotting
imagesc_gray(dataset.get_target(2),1,'Original','131');
imagesc_gray(noisy,1,'Salt & Pepper Noise','132');
imagesc_gray(denoised,1,'TV-l1 Image Denoising','133');