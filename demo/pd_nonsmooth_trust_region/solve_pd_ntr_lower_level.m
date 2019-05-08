function [y] = solve_pd_ntr_lower_level(lambda,noisy)
    %% Solving the Lower Level Problem
    param_solver.verbose = 0;
    param_solver.maxiter = 2000;
    param_solver.tol = 1e-2;

    %% Define the cell matrices
    [M,N] = size(noisy);
    id_op = IdentityOperator([M,N]);
    z = noisy;
    gradient_op = FinDiffOperator([M,N],'fn');
    q = zeros(M,N,2);
    alpha = ones(M,N);
    gamma = 0; % NO Huber regularization
    po = PatchOperator(size(lambda),[M,N]);
    lambda_out = po.val(lambda);

    %% Call the solver
    [y,~] = solve_generic_l1_l2({lambda_out},{alpha},{id_op},{gradient_op},z,{q},gamma,0*noisy,param_solver);
end

