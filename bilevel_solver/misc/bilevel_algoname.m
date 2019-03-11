function algo = bilevel_algoname(name)

  switch lower(name)
  case 'nonsmooth_trust_region'
    algo = nonsmooth_trust_region_alg();
  case 'bfgs'
    algo = bilevel_bfgs_alg();
  end
end
