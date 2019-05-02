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

id_op = IdentityOperator([M,N]);
grad_op = FinDiffOperator([M,N],'fn');
q = zeros(M,N,2);

alpha = 0.5*triu(ones(M,N))+10.5*tril(ones(M,N));
lambda = ones(M,N);

%% Solving
[denoised,gap] = solve_generic_l1_l2({lambda},{alpha},{id_op},{grad_op},noisy,q,0,0*noisy,param);