function [stop, crit, s, iter, info] = convergence_test(sol, s, iter, fg, Fp, info, param)
%CONVERGENCE_TEST Test the convergence of an algorithm

  stop = 0;
  crit = 'NOT_DEFINED';

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
  case 'rel_norm_primal'
    rel_norm_primal = norm(sol(:)-s.prev_sol(:))/norm(sol(:));
    info.rel_norm_primal(iter) = rel_norm_primal;
    if (abs(rel_norm_primal) < param.tol) && iter > 1
      crit = 'REL_NORM_PRIMAL';
      stop = 1;
    end
    s.prev_sol = sol;
  case 'rel_norm_dual'
    [rel_norm_dual,s.dual_var_old] = eval_dual_var(s.dual_var, s.dual_var_old);
    info.rel_norm_dual(iter) = rel_norm_dual;
    if (abs(rel_norm_dual) < param.tol) && iter > 1
      crit = 'REL_NORM_DUAL';
      stop = 1;
    end
  case 'rel_norm_primal_dual'
    rel_norm_primal = norm(sol(:)-s.prev_sol(:))/norm(sol(:));
    info.rel_norm_primal(iter) = rel_norm_primal;
    [rel_norm_dual, s.dual_var_old] = eval_dual_var(s.dual_var, s.dual_var_old);
    info.rel_norm_dual(iter) = rel_norm_dual;
    if (abs(rel_norm_primal) < param.tol) && (abs(rel_norm_dual) < param.tol) && iter > 1
      crit = 'REL_NORM_PRIMAL_DUAL';
      stop = 1;
    end
    s.prev_sol = sol;
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
    error('Unknown stopping criterion!')
  end

  % Handle verbosity
  if param.verbose >= 2
    switch lower(param.stopping_criterion)
    case 'rel_norm_obj'
      fprintf('  f(x^*) = %e, rel_eval = %e\n', curr_eval, rel_eval);
    case 'rel_norm_primal'
      fprintf('   Relative norm of the primal variables: %e\n', rel_norm_primal);
    case 'rel_norm_dual'
      fprintf('   Relative norm of the dual variables: %e\n', rel_norm_dual);
    case 'rel_norm_primal_dual'
      fprintf('   Relative norm of the primal variables: %e\n', rel_norm_primal);
      fprintf('   Relative norm of the dual variables: %e\n', rel_norm_dual);
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

function [rel_norm, x_old] = eval_dual_var(x1,x2)
  if nargin < 2
    if iscell(x1)
      x2 = cell(length(x1),1);
      for ii = 1:length(x1)
        x2{ii} = 0;
      end
    else
      x2 = 0;
    end
  end

  if iscell(x1)
    rel_norm = 0;
    for ii = 1:length(x1)
      tmp = norm(x1{ii}(:)-x2{ii}(:))/(norm(x1{ii}(:))+eps);
      if tmp > rel_norm
        rel_norm = tmp;
      end
    end
  else
    rel_norm = norm(x1(:)-x2(:))/(norm(x1(:))+eps);
  end
  x_old = x1;
end
