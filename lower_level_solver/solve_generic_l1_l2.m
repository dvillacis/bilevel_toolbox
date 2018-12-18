function [sol,gap] = solve_generic_l1_l2(lambda,alpha,Ks,Bs,z,q,gamma,xinit,param)
% SOLVE_GENERIC_L1_L2 Genric solver for several image processing problems
% This solver receives an abstract initial structure to support different
% image processing models and solves those by using a Chabolle-Pock algorithm.
% INPUTS
%   lambda: l2 fidelity weights
%   alpha: l1 fidelity weights
%   Ks: cell array of matrices corresponding to l2 terms
%   Bs: cell array of matrices corresponding to the l1 terms
%   z: l2 data
%   q: l1 data
%   gamma: Huber regularization parameter for l1 terms
%   param: struct with the algorithm specific parameters
% OUTPUTS
%   sol: minimizer for the optimization problem
%   gap: primal-dual gap values per iteration
%
  global FUBAR; %% TODO: Remove and rename FUBAR to PSEUDO_GAP_STATE
  FUBAR=[];

  % Start the counter
  t1 = tic;

  % Test maxiter parameter
  if ~isfield(param,'maxiter')
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

  % Test for cell input Ks
  if ~iscell(Ks)
    error('Ks must be a cell array.');
  end

  % Test for cell input Bs
  if ~iscell(Bs)
    error('Bs must be a cell array.');
  end
  
  % Test for cell input alpha
  if ~iscell(alpha)
      error('alpha must be a cell array.');
  end

  % Test for input vector
  if ~isvector(xinit)
    error('xinit must ve a vector, not a matrix.');
  end

  L = sqrt(8+1);
  tau = 0.05/L;
  sigma = 0.99/(tau*L^2);

  sol = xinit;
  sol_=sol;

  % Concatenate l2 and l1 matrices
  Kbb = cat(1,Ks{:},Bs{:});
  y = zeros(size(Kbb,1),1);

  gap = [compute_generic_l1_l2_pd_gap(sol,y,Ks,Bs,lambda,alpha,z,q)];

  if param.verbose > 1
    fprintf('generic_l1_l2: iter = %4d, gap = %f\n', 0, gap(1));
  end

  for k = 1:param.maxiter

    % Dual update
    y = y + sigma*Kbb*sol_;
    y = calc_prox(y,z,q,Ks,Bs,lambda,alpha,sigma);

    % Primal update
    sol_ = sol;
    sol = sol - tau * Kbb' * y;

    % Interpolation step
    sol_ = 2*sol - sol_;

    ga = compute_generic_l1_l2_pd_gap(sol,y,Ks,Bs,lambda,alpha,z,q);
    gap = [gap, ga];

    if mod(k, param.check) == 0 && param.verbose > 1
      fprintf('generic_l1_l2: iter = %4d, pseudo-gap = %f\n', k, ga);
    end

    % Stopping criteria
    if ga < param.tol
        fprintf('generic_l1_l2: gap-drop= %f, stopping criteria met.\n',ga);
        break;
    end

  end

  % Print summary
  if param.verbose>0
    fprintf(['\n ','GENERIC_L1_L2_CHAMBOLLE_POCK',':\n']);
    fprintf(' %i iterations\n', k);
    fprintf(' Primal-Dual Pseudo-Gap: %f \n', gap(end));
    fprintf(' Execution Time: %f \n\n', toc(t1));
  end

end

function prox = calc_prox(y,z,q,Ks,Bs,lambda,alpha,tau)
    index = 0;
    for k=1:length(Ks)
        n = size(Ks{k},1);
        y(index+1:index+n) = (y(index+1:index+n)-tau.*z)./(1+0.5*tau*(1./lambda)); % l2 proximal
        index = index + n;
    end
    for l = 1:length(Bs)
        n = size(Bs{l},1);
        y(index+1:index+n) = projection_l2_ball(y(index+1:index+n)-tau*q,alpha{l});
        index = index + n;
    end
    prox = y;
end
