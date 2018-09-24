function algo = bilevel_algoname(name)

  switch lower(name)
  case 'nonsmooth_trust_region'
    algo = nonsmooth_trust_region_alg();
  case 'semismooth_newton'
    algo = semismooth_newton_alg();
  end
end
