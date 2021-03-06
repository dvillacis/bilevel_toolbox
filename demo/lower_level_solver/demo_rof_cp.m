clear all;
clc;

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset','*_circle_original.png','*_circle_noisy.png');

%% Load image
noisy = dataset.get_corrupt(2);

%% Solver Parameters
param.verbose = 2;
param.maxiter = 5000;
param.check = 200;
param.alpha = 60.0;

%% Solving
[denoised,gap] = solve_rof_cp_single_gaussian(noisy,param);

%% Plotting
figure
imagesc_gray(denoised);