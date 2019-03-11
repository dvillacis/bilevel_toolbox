function s = bilevel_bfgs_alg()
  s.name = 'BFGS';
  s.initialize = @(initial_parameter, lower_level_problem, upper_level_problem, param) bilevel_bfgs_initialize(initial_parameter,lower_level_problem,upper_level_problem,param);
  s.algorithm = @(initial_parameter,lower_level_problem,upper_level_problem,sol,s,param) bilevel_bfgs_algorithm(lower_level_problem,upper_level_problem,sol,s,param);
  s.finalize = @(initial_parameter,lower_level_problem,upper_level_problem,sol,s,param) sol;
end

function [sol,state,param] = bilevel_bfgs_initialize(initial_parameter,lower_level_problem,upper_level_problem,param)
  sol = initial_parameter;
  state.sol_history = sol;
  state.res = 1;
  state.u_history = lower_level_problem.solve(sol,upper_level_problem.dataset);
  state.BFGS = eye(size(state.u_history(1),1)); %Initial BFGS Matrix

  % Default param values
  if ~isfield(param, 'armijo_c')
      param.armijo_c = 1e-4;
  end

  if ~isfield(param, 'wolfe_c')
      param.wolfe_c = 0.1;
  end

end

function [sol,state] = bilevel_bfgs_algorithm(lower_level_problem,upper_level_problem,sol,state,param)

  % Solving the state equation (lower level solver)
  u = lower_level_problem.solve(sol,upper_level_problem.dataset);
  state.u_history = cat(3,state.u_history,u);
  % Getting current cost
  cost = upper_level_problem.eval(u,sol,upper_level_problem.dataset);

  % Saving cost history
  if ~isfield(s, 'l2_cost_history')
      state.l2_cost_history = cost;
  else
      state.l2_cost_history = [s.l2_cost_history cost];
  end

  while state.res >= param.tol
    state.res = 0;
  end

end
