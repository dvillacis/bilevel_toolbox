function stop_crit = select_stopping_criterion(algo)

  if ischar(algo)
    switch lower(algo)
      case 'admm'
        stop_crit = 'rel_norm_primal_dual';
      case 'chambolle_pock'
        stop_crit = 'rel_norm_primal_dual';
      case 'sdmm'
        stop_crit = 'rel_norm_primal_dual';
      case 'pocs'
        stop_crit = 'obj_threshold';
      case 'fb_based_primal_dual'
        stop_crit = 'rel_norm_primal_dual';
      case 'fbf_primal_dual'
        stop_crit = 'rel_norm_primal_dual';
      otherwise
        stop_crit = 'rel_norm_obj';
    end
  else
    stop_crit = 'rel_norm_obj';
  end
end
