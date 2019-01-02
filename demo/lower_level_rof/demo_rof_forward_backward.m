clear all;
close all;
clc;

init_bilevel_toolbox();

%% Load dataset
dataset = DatasetInFolder('data/playing_cards','*_playing_cards_original.tif','*_playing_cards_noisy.tif');

%% Load input image
original = dataset.get_target(5);
noisy = dataset.get_corrupt(5);

%% Solving the Lower Level Problem
param_solver.verbose = 2;
param_solver.maxiter = 900;
param_solver.alpha = 0.1;

[sol,gap] = solve_rof_fb_single_gaussian(noisy,param_solver);

%% Plotting the solution
figure(1)
subplot(1,3,1)
imagesc_gray(original,1,'Original Image');
subplot(1,3,2)
imagesc_gray(noisy,1,'Noisy Image');
subplot(1,3,3)
imagesc_gray(sol,1,'Denoised Image');

figure(2)
loglog(gap);
