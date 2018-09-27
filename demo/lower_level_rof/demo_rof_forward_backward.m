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

%imagesc_gray(original,1,'Original Image');
%imagesc_gray(noisy,2,'Noisy Image');

%% Setting Lower Level Solver
f2.grad = @(u) u - noisy;
f2.eval = @(u) 0.5 * norm(u-noisy)^2;
f2.beta = 1;

param_tv.verbose = verbose - 1;
alpha = 0.1;
f1.prox = @(u,T) prox_tv(u,alpha*T,param_tv);
f1.eval = @(u) alpha*norm_tv(u);

%% Solving the Lower Level Problem
param_solver.verbose = verbose;
param_solver.maxit = 300;
param_solver.tol = 1e-7;
param_solver.gamma = 0.1;

u = zeros(size(noisy));

[sol,info] = solvep(u,{f1,f2},param_solver);

%% Plotting the solution
imagesc_gray(original,1,'Original Image');
imagesc_gray(noisy,2,'Noisy Image');
imagesc_gray(sol,3,'Denoised Image');
