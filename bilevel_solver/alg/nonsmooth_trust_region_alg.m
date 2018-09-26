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
  [lower_sol, ~] = lower_level_problem.solve(noisy,sol);

  % Solving the adjoint state
  [p,grad] = adjoint_solver(lower_sol,original,sol);

end

function nXi = xi(p,m,n)
  p = reshape(p,m*n,2);
  a = sqrt(sum(p.^2,2));
  nXi = [a;a];
end

function prod = outer_product(p,q,m,n)
  p = reshape(p,m*n,2);
  q = reshape(q,m*n,2);
  a = p(:,1).*q(:,1);
  b = p(:,1).*q(:,2);
  c = p(:,2).*q(:,1);
  d = p(:,2).*q(:,2);
  prod = [spdiags(a,0,m*n,m*n) spdiags(b,0,m*n,m*n); spdiags(c,0,m*n,m*n) spdiags(d,0,m*n,m*n)];
end

function [adj,grad] = adjoint_solver(u,original,sol)
  % Get the adjoint state
  [m,n] = size(u);
  nabla = gradient_matrix(m,n);
  Ku = nabla*u(:);
  nKu = xi(Ku,m,n);
  act = (nKu<1e-7);
  inact = 1-act;
  Act = spdiags(act,0,2*m*n,2*m*n);
  Inact = spdiags(inact,0,2*m*n,2*m*n);
  denominador = Inact*nKu+act;
  prodKuKu = outer_product(Ku./(denominador.^3),Ku,m,n);
  A = speye(m*n);
  B = nabla';
  C = sol*Inact*(prodKuKu-spdiags(1./denominador,0,2*m*n,2*m*n))*nabla;
  D = speye(2*m*n);
  E = Act*nabla;
  F = sparse(2*m*n,2*m*n);
  Adj = [A B;C D;E F];
  Track = [u(:)-original(:);sparse(4*m*n,1)];
  mult = Adj\Track;
  adj = mult(1:n*m);
  % Calculating the gradient
  Kp = nabla*adj;
  aux = Inact*(Ku./denominador);
  grad = sol*0.001*m*n - aux'*Kp;
end
