function [sol,gap] = solve_rof_fb_single_gaussian(f,param)

  % Start the counter
  t1 = tic;

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
    param.tol = 1e-4;
  end


  [M, N] = size(f);
  f = f(:);

  % Generate gradient matrix for this image
  gradient_op = FinDiffOperator([M,N])
  nabla = gradient_matrix(M,N);

  p = zeros(M*N*2,1);
  p_old = p;

  L = 8;
  gap = [];

  for k = 1:param.maxiter

    p_old = p;
    p = p - 0.99/L*(nabla*(nabla'*p - f)); %% TODO: Evitar division (precalcular afuera)
    p = reshape(p,M*N,2);
    p = reshape(bsxfun(@rdivide,p,max(1, sqrt(sum(p.^2,2))/param.alpha)), M*N*2,1);
    grad = p_old-p;
    p = p_old - 1.9*grad;

    sol = f - nabla'*p;
    [ga, ~, ~] = compute_rof_pd_gap(nabla, sol, p, f, param.alpha, 0, M, N);

    gap = [gap, ga];

    if mod(k, param.check) == 0 && param.verbose > 1
      fprintf('rof_fb: iter = %4d, gap = %f\n', k, ga);
    end

    % if ga <= param.tol
    %   break;
    % end

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
