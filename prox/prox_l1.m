function [sol,info] = prox_l1(x, gamma, param)
%PROX_L1 Proximal operator with L1 norm

  % Start the time counter
  t1 = tic;

  % Optional input arguments
  if nargin<3, param = struct; end

  % Optional input arguments
  if ~isfield(param, 'verbose'), param.verbose = 1; end
  if ~isfield(param, 'nu'), param.nu = 1; end
  if ~isfield(param, 'tol'), param.tol = 1e-3; end
  if ~isfield(param, 'maxit'), param.maxit = 200; end
  if ~isfield(param, 'At'), param.At = @(x) x; end
  if ~isfield(param, 'A'), param.A = @(x) x; end
  if ~isfield(param, 'weights'), param.weights = 1; end
  if ~isfield(param, 'pos'), param.pos = 0; end

  % Test the parameters
  if test_gamma(gamma)
    sol = x;
    info.alfo = mfilename;
    info.iter = 0;
    info.final_eval = 0;
    info.crit = '--';
    info.time = toc(t1);
  end

  param.weights = test_weights(param.weights);
  temp = param.A(x);
  if ~isfield(param, 'y'), param.y = zeros(size(temp)); end
  temp2 = param.y + soft_threshold(temp -  param.y , gamma*param.nu*param.weights) - temp;
  sol = x + 1/param.nu * param.At(temp2);
  crit = 'REL_OBJ'; iter = 1;
  dummy = temp2 +temp;
  norm_l1 = gamma*sum(param.weights(:).*abs(dummy(:)));

  % Log after the prox l1
  if param.verbose >= 1
    fprintf(['  prox_L1: ||A x-y||_1 = %e,', ' %s, iter = %i\n'], norm_l1, crit, iter);
  end

  info.algo = mfilename;
  info.iter = iter;
  info.final_eval = norm_l1;
  info.crit = crit;
  info.time = toc(t1);

end
