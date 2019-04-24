function [y] = solve_regularized_bfgs_lower_level(lambda,noisy)
    %% Solving the Lower Level Problem
    param_solver.verbose = 0;
    param_solver.maxiter = 2000;
    param_solver.tol = 1e-3;

    %% Define the cell matrices
    [M,N] = size(noisy);
    id_op = IdentityOperator([M,N]);
    z = noisy;
    gradient_op = FinDiffOperator([M,N],'fn');
    q = zeros(M,N,2);
    alpha = ones(M,N);
    gamma = 0; % NO Huber regularization

    %% Call the solver
    [y,~] = solve_generic_l1_l2({lambda},{alpha},{id_op},{gradient_op},z,q,gamma,0*noisy,param_solver);
end

