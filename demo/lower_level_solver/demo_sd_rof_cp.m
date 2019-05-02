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
param.check = 500;
param.alpha = 0.1*triu(ones(M,N))+10.5*tril(ones(M,N));

%% Solving
[denoised,gap] = solve_sd_rof_cp_single_gaussian(noisy,param);