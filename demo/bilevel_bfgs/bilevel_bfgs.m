clear all;
close all;
clc;

init_bilevel_toolbox();

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset_single_gaussian','*_circle_original.png','*_circle_noisy.png');

% Define lower level problem
lower_level_problem.solve = @(lambda,dataset) solve_lower_level(lambda,dataset.get_corrupt(1));

% Define upper level problem
upper_level_problem.eval = @(u,lambda,dataset) eval_upper_level(u,lambda,dataset);
upper_level_problem.gradient = @(u,lambda,dataset,params) solve_gradient(u,lambda,dataset,params);
upper_level_problem.dataset = dataset;

%% Solving the bilevel problem
bilevel_param.verbose = 2;
bilevel_param.maxit = 50;
bilevel_param.tol = 1e-3;
bilevel_param.algo = 'NONSMOOTH_TRUST_REGION';
bilevel_param.armijo_c = 1e-4;
bilevel_param.wolfe_c = 0.1;
lambda = 3;
[sol,info] = solve_bilevel(lambda,lower_level_problem,upper_level_problem,bilevel_param);


%% Auxiliary functions
function cost = eval_upper_level(u,lambda,dataset)
    original = dataset.get_target(1);
    cost = 0.5*norm(u(:)-original(:)).^2 + 0.0001*norm(lambda);
end

function u = solve_lower_level(lambda,noisy)
    param_solver.verbose = 0;
    param_solver.maxiter = 2000;
    param_solver.tol = 1e-2;
    %% Define the cell matrices
    [M,N] = size(noisy);
    id_op = IdentityOperator([M,N]);
    z = noisy;
    gradient_op = FinDiffOperator([M,N],'fn');
    B = gradient_op.matrix();
    q = zeros(M,N,2);
    alpha = 1;
    gamma = 0;
    [u,~] = solve_generic_l1_l2({lambda},{alpha},{id_op},{gradient_op},z,q,gamma,0*noisy,param_solver);
end
