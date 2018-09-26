function s = nonsmooth_trust_region_alg()
  s.name = 'NONSMOOTH_TRUST_REGION';
  s.initialize = @(x_0, dataset, lower_level_problem, upper_level_problem, param) nonsmooth_trust_region_initialize(x_0,dataset,lower_level_problem,upper_level_problem,param);
  s.algorithm = @(x_0,dataset,lower_level_problem,upper_level_problem,sol,s,param) nonsmooth_trust_region_algorithm(dataset,lower_level_problem,upper_level_problem,sol,s,param);
  s.finalize = @(x_0,lower_level_problem,upper_level_problem,sol,s,param) sol;
end

function [sol,s,param] = nonsmooth_trust_region_initialize(x_0,dataset,lower_level_problem,upper_level_problem,param)

  s.x_n = {};
  sol = x_0;

  % Test if lower level problem has a solve method
  if ~isfield(lower_level_problem, 'solve')
    error('Lower Level Problem struct does not provide a solve method.')
  end

  % TODO: Verify if dataset is correct

end

function [sol,s] = nonsmooth_trust_region_algorithm(dataset,lower_level_problem,upper_level_problem,sol,s,param)

  % Load dataset
  original = dataset.get_target(1);
  noisy = dataset.get_corrupt(1);

  % Solving the state equation (lower level solver)
  [lower_sol, lower_info] = lower_level_problem.solve(noisy,sol);

  % Get the adjoint state
  [Kux,Kuy] = gradient_op(lower_sol);
  Ku = cat(3,Kux,Kuy);
  nKu = sqrt(Kux.^2+Kuy.^2); % Pixelwise l2 norm of the dual variable Ku
  Act = (nKu<1e-2);
  Inact = 1-Act;
  denominador = Inact.*nKu+Act;
  prodKuKu = outer_product(Ku./(denominador.^3),Ku);

end
