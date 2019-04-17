function [sol,gap] = solve_rof_fb_single_gaussian(f,param)

  % Test alpha parameter
  if ~isfield(param,'alpha')
    error('Lower Level Problem struct does not provide an ALPHA parameter.')
  end

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
    param.tol = 1e-3;
  end


  [M, N] = size(f);

  % Generate gradient matrix for this image
  op = FinDiffOperator([M,N],'fn');

  p = zeros(M,N,2);
  p_old = p;

  L = 8;
  gap = [];

  % Auxiliary terms
  a = 0.99/L;
  b = 1/param.alpha;

  % Start the counter
  t1 = tic;

  for k = 1:param.maxiter

    p_old = p;
    p = p - a * op.val(op.conj(p)-f);
    p = bsxfun(@rdivide,p,max(1,b*rssq(p,3)));
    
    %p = p - a*(nabla*(nablat*p - f));
    %p = reshape(p,M*N,2);
    %p = reshape(bsxfun(@rdivide,p,max(1, b*rssq(p,2))), M*N*2,1);
    
    grad = p_old-p;
    p = p_old - 1.9*grad;

    sol = f - op.conj(p);
    [ga, ~, ~] = compute_rof_pd_gap(op.val(sol), sol, op.conj(p), f, param.alpha);

    gap = [gap, ga];

    if mod(k, param.check) == 0 && param.verbose > 1
      fprintf('rof_fb: iter = %4d, gap = %f\n', k, ga);
    end

    if ga <= param.tol
      break;
    end

  end
  sol = reshape(sol,M,N);

  % Print summary
  if param.verbose>0
    fprintf(['\n ','ROF_FORWARD_BACKWARD',':\n']);
    fprintf(' %i iterations\n', k);
    fprintf(' Primal-Dual Gap: %f \n\n', gap(end));
    fprintf(' Execution Time: %f \n\n', toc(t1));
  end

end

function [gap, primal, dual] = compute_rof_pd_gap(val_u, u, conj_p, f, alpha)

  norm_val_u = rssq(val_u,3);
  
  primal = alpha * sum(norm_val_u(:)) + 0.5 * sumsqr(u(:)-f(:));
  dual   = f(:)'*conj_p(:) - 0.5 * sumsqr(conj_p(:));
  
  gap = primal-dual;
end
