function [info, iter, s] = bilevel_initialize_convergence_variable(sol, s, lower_level_problem, upper_level_problem, param)

  switch lower(param.stopping_criterion)
    case 'rel_norm_obj'
      info.objective = nan(param.maxit+1,1);
      info.rel_eval = nan(param.maxit,1);
      info.objective(1) = eval_function(fg, Fp, sol, s, param);
    case 'obj_increase'
      info.objective = nan(param.maxit+1,1);
      info.objective(1) = eval_function(fg, Fp, sol, s, param);
    case 'obj_threshold'
      info.objective = nan(param.maxit+1,1);
      info.objective(1) = eval_function(fg, Fp, sol, s, param);
    otherwise
      error('Unknown stopping criterion!');
  end
end
