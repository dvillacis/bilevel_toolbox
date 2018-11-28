function [sol,gap] = solve_generic_l1_l2(alpha,beta,As,Bs,omega,zeta,gamma,f,param)
% SOLVE_GENERIC_L1_L2 Genric solver for several image processing problems
% This solver receives an abstract initial structure to support different
% image processing models and solves those by using a Chabolle-Pock algorithm.
% INPUTS
%   alpha: l2 fidelity weights
%   beta: l1 fidelity weights
%   As: cell array of matrices corresponding to l2 problem
%   Bs: cell array of matrices corresponding to the l1 problem
%   gamma: Huber regularization parameter
%   f: contaminated image
%   param: struct with the algorithm specific parameters
% OUTPUTS
%   sol: minimizer for the optimization problem
%   gap: primal-dual gap values per iteration
%

  % Start the counter
  t1 = tic;

  % Test maxiter parameter
  if ~isfield(param,'maxit')
    param.maxiter = 1000;
  end

  % Test check parameter
  if ~isfield(param,'check')
    param.check = 100;
  end

  % Test verbose parameter
  if ~isfield(param,'verbose')
    param.verbose = 1;
  end

  % Test tol parameter
  if ~isfield(param,'tol')
    param.tol = 1e-4;
  end



  % Print summary
  if param.verbose>0
    fprintf(['\n ','GENERIC_L1_L2_CHAMBOLLE_POCK',':\n']);
    fprintf(' %i iterations\n', k);
    fprintf(' Primal-Dual Gap: %f \n\n', gap(end));
    fprintf(' Execution Time: %f \n\n', toc(t1));
  end

end
