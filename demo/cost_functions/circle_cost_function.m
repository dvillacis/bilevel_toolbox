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
lower_level_problem.solve = @(alpha) solve_lower_level(alpha,noisy);

c = [];
r = 0.01:0.001:0.2;
i=1;
for a = r
    sol = lower_level_problem.solve(a);
    cost = 0.5*norm(original(:)-sol(:)).^2;
    c = [c cost];
    if mod(i,31)==0
        fprintf('Finished %.3f with cost %f\n',a,cost);
    end
    i=i+1;
end

plot(r,c)

function y = solve_lower_level(alpha,noisy)
  param_lower_level.maxit = 1000;
  param_lower_level.alpha = alpha;
  param_lower_level.verbose = 0;
  y = solve_rof_cp_single_gaussian(noisy,param_lower_level);
end