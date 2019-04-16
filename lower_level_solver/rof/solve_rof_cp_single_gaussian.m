function [sol,gap] = solve_rof_cp_single_gaussian(f,param)

  % Test alpha parameter
  if ~isfield(param,'alpha')
    error('Lower Level Problem struct does not provide an ALPHA parameter.')
  end

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


  [M, N] = size(f);
  %f = f(:);

  op = FinDiffOperator([M,N],'fn');

  p = zeros(M,N,2);
  sol = f;
  sol_ = sol;
  L = sqrt(8);
  tau = 0.01;
  sigma = 1/tau/L^2;
  gap = [];

  % Auxiliary terms
  a = 1/(1+tau);
  b = 1/param.alpha;

  % Start the counter
  t1 = tic;

  for k = 1:param.maxiter

    val_u = op.val(sol_);
    p = p + sigma*val_u;
    p = bsxfun(@rdivide,p,max(1,b*rssq(p,3)));
    %p = reshape(p,M*N,2);
    %p = reshape(bsxfun(@rdivide,p,max(1, b*rssq(p,2))), M,N);

    conj_p = op.conj(p);
    sol_ = sol;
    sol = sol - tau*conj_p;
    sol = a*(sol+tau*f);
    sol_ = 2*sol -sol_;

    ga = compute_rof_pd_gap(val_u, sol, conj_p, f, param.alpha);

    gap = [gap, ga];

    if mod(k, param.check) == 0 && param.verbose > 1
      fprintf('rof_cp: iter = %4d, gap = %f\n', k, ga);
    end

    if ga < param.tol
      break;
    end

  end
  sol = reshape(sol,M,N);

  % Print summary
  if param.verbose>0
    fprintf(['\n ','ROF_CHAMBOLLE_POCK',':\n']);
    fprintf(' %i iterations\n', k);
    fprintf(' Primal-Dual Gap: %f \n', gap(end));
    fprintf(' Execution Time: %f \n\n', toc(t1));
  end

end
