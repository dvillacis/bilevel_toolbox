clear all;
clc;

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset_single_gaussian','*_circle_original.png','*_circle_noisy.png');

%% Load image
noisy = dataset.get_corrupt(1);

%% Solver Parameters
param.verbose = 2;
param.maxiter = 2000;
param.check = 200;

id_op = IdentityOperator(size(noisy));
grad_op = FinDiffOperator(size(noisy),'fn');
q = zeros(size(noisy,1),size(noisy,2),2);

alpha = 1;
lambda = 4.5;

%% Solving
[denoised,gap] = solve_generic_l1_l2({lambda},{alpha},{id_op},{grad_op},noisy,q,0,0*noisy,param);
