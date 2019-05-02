clear all;
clc;

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset','*_circle_original.png','*_circle_noisy.png');

%% Load image
noisy = dataset.get_corrupt(1);

%% Solver Parameters
param.verbose = 2;
param.maxiter = 2000;
param.check = 200;
param.alpha = 0.1;

%% Solving
[denoised,gap] = solve_rof_fb_single_gaussian(noisy,param);