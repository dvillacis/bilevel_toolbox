function [sol,info] = solve_bilevel(x_0, lower_level_problem, upper_level_problem, param)

  % Start the counter
  t1 = tic;

  % Check for optional input arguments
  if nargin<4, param=struct; end

  % Check mandaotry arguments
  if nargin<3
    error('Not enough input parameters');
  end

  % Default param values
  if ~isfield(param, 'tol'), param.tol=1e-4 ; end
  if ~isfield(param, 'maxit'), param.maxit=200; end
  if ~isfield(param, 'verbose'), param.verbose=1 ; end

  % TODO: test for lower level correctness


  % Select the stopping criterion
  if ~isfield(param, 'stopping_criterion')
    param.stopping_criterion = bilevel_select_stopping_criterion(param.algo);
  end

  % TODO: Setup algorithm selection strategy
  algo = get_bilevel_algo(param.algo);

  % Initialization
  [sol,s,param] = algo.initialize(x_0,lower_level_problem,upper_level_problem,param);
  [info,iter,s] = bilevel_initialize_convergence_variable(sol,s,lower_level_problem,upper_level_problem,param);

  % Main Loop
  while 1

    if param.verbose >= 1
      fprintf('Bilevel Iter %.3i: ',iter);
    end

    [sol,s] = algo.algorithm(x_0,lower_level_problem,upper_level_problem,sol,s,param);

    [stop,crit,s,iter,info] = bilevel_convergence_test(sol,s,iter,lower_level_problem,upper_level_problem,info,param);

    % [sol,param] = bilevel_post_process(sol,iter,info,param);

    if stop, break; end

  end

  sol = algo.finalize(x_0,lower_level_problem,upper_level_problem,sol,s,param);

  info.algo = algo.name;
  info.iter = iter;
  info.crit = crit;
  info.time = toc(t1);
  
  % Print summary
  if param.verbose>0
    fprintf(['\n ','%s',':\n'],info.algo);
    fprintf(' %i iterations\n', info.iter);
    %fprintf(' Primal-Dual Gap: %f \n', gap(end));
    fprintf(' Execution Time: %f \n\n', info.time);
  end

end

function algo = get_bilevel_algo(name)
  if isstruct(name)
    algo = name;
  elseif ischar(name)
    algo = bilevel_algoname(name);
  else
    error('The algorithm is not a string nor a struct');
  end
end
