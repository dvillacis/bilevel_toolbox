function [sol,info] = solvep(x_0, F, param)

  % Start the counter
  t1 = tic;

  % Check for optional input arguments
  if nargin<3, param=struct; end

  % Check mandatory arguments
  if nargin<2
    error('Not enough input arguments')
  end

  % Default param values
  if ~isfield(param, 'tol'), param.tol=1e-4 ; end
  if ~isfield(param, 'maxit'), param.maxit=200; end
  if ~isfield(param, 'verbose'), param.verbose=1 ; end
  if ~isfield(param, 'lambda'), param.lambda=0.99 ; end
  if ~isfield(param, 'fast_eval'), param.fast_eval = 0  ; end
  if ~isfield(param, 'debug_mode'), param.debug_mode = 0 ; end

  % Test input for eval method
  if ~iscell(F)
    F = {F};
  end
  F = test_eval(F);

  % Test input for grad and prox
  [Fg,Fp] = prepare_function(F,param);
  if ~isfield(param,'gamma')
    if numel(Fg)
      param.gamma = compute_gamma(Fg);
    else
      param.gamma = 1;
    end
  else
    if param.verbose >= 1
      fprintf('The time step is set manually to %g\n', param.gamma);
    end
  end

  % Select the algorithm
  if ~isfield(param, 'algo'), param.algo = select_solver(Fg,Fp); end

  % Select the stopping criterion
  if ~isfield(param, 'stopping_criterion')
    param.stopping_criterion = select_stopping_criterion(param.algo);
  end

  algo = get_algo(param.algo);

  % Transform all smooth functions into one function
  fg = add_smooth_functions(Fg);

  if param.verbose >= 1
    fprintf(['Algorithm selected: ', algo.name,' \n']);
  end

  [sol,s,param] = algo.initialize(x_0, fg, Fp, param);

  % Initialization
  [info,iter,s] = initialize_convergence_variable(sol, s, fg, Fp, param);

  % Main Loop
  while 1

    if para.verbose >= 1
      fprintf('Iter %.3i: ', iter);
    end

    [sol,s] = algo.algorithm(x_0,fg,Fp,sol,s,param);

    [stop,crit,s,iter,info] = convergence_test(sol,s,iter,fg,Fp,info,param);

    [sol,param] = post_process(sol,iter,info,param);

    if stop, break; end

  end

  sol = algo.finalize(x_0,fg,Fp,sol,s,param);

  summary_print(s,info,iter,algo,crit,param);

  info.algo = algo.name;
  info.iter = iter;
  info.crit = crit;
  info.time = toc(t1);

  % Return the dual variable if available
  if isfield(s,'dual_var')
    info.dual_var = s.dual_var;
  end

end

function gamma = compute_gamma(Fg)
  beta = 0;
  for ii = 1:length(Fg)
    beta = beta + Fg{ii}.beta;
  end
  if beta == 0
    gamma = 1;
  elseif beta > 0
    gamma = 1/beta;
  else
    error('Error in the lipschitz constant');
  end
end

function algo = get_algo(name)
  if isstruct(name)
    algo = name;
  elseif ischar(name)
    algo = algoname(name);
  else
    error('The algorithm is not a string nor a struct!');
  end
end
