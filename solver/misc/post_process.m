function [sol,param] = post_process(sol, iter, info, param)

  if ~isfield(param,'gamma'), param.gamma = 1; end
  if ~isfield(param,'do_sol'), param.do_sol = @(x) x.sol; end
  if ~isfield(param,'do_ts'), param.do_ts = @(x) x.gamma; end

  info_iter.sol = sol;
  info_iter.iter = iter;
  info_iter.gamma = param.gamma;
  info_iter.info = info;

  sol = param.do_sol(info_iter);
  param.gamma = param.do_ts(info_iter);

end
