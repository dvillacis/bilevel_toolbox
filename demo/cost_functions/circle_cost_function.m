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
r = 5:0.5:11;
i=1;
for a = r
    sol = lower_level_problem.solve(a);
    cost = 0.5*norm(original(:)-sol(:)).^2;
    c = [c cost];
    if mod(i,3)==0
        fprintf('Finished %.3f with cost %f\n',a,cost);
    end
    i=i+1;
end

plot(r,c)

function sol = solve_lower_level(lambda,noisy)
    %% Solving the Lower Level Problem
    param_solver.verbose = 0;
    param_solver.tol = 1e-5;

    %% Define the cell matrices
    [M,N] = size(noisy);
    K = speye(M*N);
    z = noisy(:);
    %lambda = 3;
    gradient = FinDiffOperator([M,N],'fn');
    B = gradient.matrix();
    q = zeros(2*M*N,1);
    alpha = 1;

    gamma = 0; % NO Huber regularization

    %% Call the solver
    [sol,~] = solve_generic_l1_l2(lambda,{alpha},{K},{B},z,q,gamma,0*noisy(:),param_solver);
end