function stop_crit = bilevel_select_stopping_criterion(algo)

  if ischar(algo)
    switch lower(algo)
      case 'nonsmooth_trust_region'
        stop_crit = 'rel_norm_obj';
      otherwise
        stop_crit = 'obj_threshold';
    end
  else
    stop_crit = 'rel_norm_obj';
  end
end
