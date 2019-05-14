function [u] = solve_ntr_lower_level(lambda,noisy)
    param_solver.verbose = 0;
    param_solver.maxiter = 100;
    param_solver.tol = 1e-2;
    %% Define the cell matrices
    [M,N] = size(noisy);
    id_op = IdentityOperator([M,N]);
    z = noisy;
    gradient_op = FinDiffOperator([M,N],'fn');
    q = zeros(M,N,2);
    alpha = 1;
    gamma = 0;
    [u,~] = solve_generic_l1_l2({lambda},{alpha},{id_op},{gradient_op},z,{q},gamma,0*noisy,param_solver);
end
