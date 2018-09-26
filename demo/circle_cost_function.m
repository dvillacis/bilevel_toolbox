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

%% Setup for the lower level problem
lower_level_problem.solve = @(noisy,alpha) solve_rof_single_gaussian_lower_level(noisy,alpha);

c = [];
r = 0.01:0.01:0.15;

for a = r
    [sol,info] = lower_level_problem.solve(noisy,a);
    cost = 0.5*norm(original(:)-sol(:)).^2;
    c = [c cost];
    fprintf('Finished %.3f with cost %f in %d iterations\n',a,cost,info.iter);
end

plot(r,c)