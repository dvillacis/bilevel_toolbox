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
  f = f(:);

  gradient = FinDiffOperator([M,N],'fn');
  nabla = gradient.matrix();
  nablat = nabla';
  %nabla = gradient_matrix(M,N);

  p = zeros(M*N*2,1);
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

    p = p + sigma*nabla*sol_;
    p = reshape(p,M*N,2);
    p = reshape(bsxfun(@rdivide,p,max(1, b*rssq(p,2))), M*N*2,1);

    sol_ = sol;
    sol = sol - tau*nablat*p;
    sol = a*(sol+tau*f);
    sol_ = 2*sol -sol_;

    ga = compute_rof_pd_gap(nabla, nablat, sol, p, f, param.alpha, M, N);

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
