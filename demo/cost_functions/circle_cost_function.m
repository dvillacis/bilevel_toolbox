clear all;
close all;
clc;

verbose = 2;

init_bilevel_toolbox();

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset_single_gaussian','*_circle_original.png','*_circle_noisy.png');
%dataset = DatasetInFolder('data/cameraman_dataset_single_gaussian','*_cameraman_original.png','*_cameraman_noisy.png');

%% Load input image
original = dataset.get_target(1);
noisy = dataset.get_corrupt(1);

%% Setup for the lower level problem
lower_level_problem.solve = @(lambda) solve_lower_level(lambda,noisy);

r = 9.5:0.01:9.8;
l2_vals = zeros(length(r),1);
ssim_vals = zeros(length(r),1);
psnr_vals = zeros(length(r),1);

i=1;
for a = r
    sol = lower_level_problem.solve(a);
    cost = 0.5*norm(original(:)-sol(:)).^2 + 0.0001*0.5*norm(a).^2;
    s = ssim_index(original,sol);
    p = psnr(original,sol,1);
    l2_vals(i) = cost;
    ssim_vals(i) = s;
    psnr_vals(i) = p;
    if mod(i,4)==0
        fprintf('Finished %.3f with cost %f\n',a,cost);
    end
    i=i+1;
end

%% Plotting
figure
plot(r,l2_vals);

figure
plot(r,ssim_vals);

figure
plot(r,psnr_vals);

%% Auxiliary functions
function sol = solve_lower_level(lambda,noisy)
    % Solving the Lower Level Problem
    param_solver.verbose = 1;
    param_solver.maxiter = 3000;

    % Define the cell matrices
    [M,N] = size(noisy);
    K = speye(M*N);
    z = noisy(:);
    %lambda = 3;
    gradient = FinDiffOperator([M,N],'fn');
    B = gradient.matrix();
    q = zeros(2*M*N,1);
    alpha = 1;

    gamma = 0; % NO Huber regularization

    % Call the solver
    [sol,~] = solve_generic_l1_l2(lambda,{alpha},{K},{B},z,q,gamma,0*noisy(:),param_solver);
    sol = reshape(sol,M,N);
end