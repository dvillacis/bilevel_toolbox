clear all;
clc;

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset_single_gaussian','*_circle_original.png','*_circle_noisy.png');

%% Load image
noisy = dataset.get_corrupt(1);
[M,N] = size(noisy);

%% Solver Parameters
param.verbose = 2;
param.maxiter = 2000;
param.check = 200;

id_op = IdentityOperator([M,N]);
grad_op = FinDiffOperator([M,N],'fn');
q = zeros(M,N,2);

alpha = 1;
lambda = 4.5;

%% Solving
[denoised,gap] = solve_generic_l1_l2({lambda},{alpha},{id_op},{grad_op},noisy,q,0,0*noisy,param);
