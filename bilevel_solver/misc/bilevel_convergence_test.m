function [stop, crit, s, iter, info] = bilevel_convergence_test(sol, s, iter, lower_level_problem, upper_level_problem)

  stop = 0;
  crit = 'NOT DEFINED';

  switch lower(param.stopping_criterion)
    case 'rel_norm_obj'
      curr_eval = eval_objective(fg,Fp,sol,s,param);
      info.objective(iter+1) = curr_eval;
      rel_eval = abs(curr_eval - info.objective(iter))/(curr_eval + eps);
      info.rel_eval(iter) = rel_eval;
      if (abs(rel_eval) < param.tol) && iter > 1
        crit = 'REL_NORM_OBJ';
        stop = 1;
      end
    case 'obj_increase'
      curr_eval = eval_objective(fg, Fp, sol, s, param);
      info.objective(iter+1) = curr_eval;
      if curr_eval >= info.objective(iter)
        crit = 'OBJ_INCREASE';
        stop = 1;
      end
    case 'obj_threshold'
      curr_eval = eval_objective(fg, Fp, sol, s, param);
      info.objective(iter+1) = curr_eval;
      if (curr_eval < param.tol) && iter > 1
        crit = 'OBJ_THRESHOLD';
        stop = 1;
      end
    otherwise
      error('Unknown stopping criterion!');
  end

  % Handle verbosity
  if param.verbose >= 2
    switch lower(param.stopping_criterion)
    case 'rel_norm_obj'
      fprintf('  f(x^*) = %e, rel_eval = %e\n', curr_eval, rel_eval);
    case 'obj_increase'
      fprintf('  f(x^*) = %e, prev_it:  %e\n', curr_eval, info.objective(iter));
    case 'obj_threshold'
      fprintf('  f(x^*) = %e\n', curr_eval);
    otherwise
      error('Unknown stopping criterion!');
    end
  end

  % Stopping if too many iterations
  if iter >= param.maxit
    crit = 'MAX_IT';
    stop = 1;
  end

  % Performed updates
  if ~stop
    iter = iter + 1;
  end

end

function curr_eval = eval_objective(fg, Fp, sol, s, param)
  curr_eval = eval_function(fg, Fp, sol, s, param);
  if ~(numel(curr_eval)==1)
    error('A least, one of your evaluation function does not return a scalar.');
  end
end
