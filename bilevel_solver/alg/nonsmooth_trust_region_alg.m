function s = nonsmooth_trust_region_alg()
  s.name = 'NONSMOOTH_TRUST_REGION';
  s.initialize = @(x_0, lower_level_problem, upper_level_problem, param) nonsmooth_trust_region_initialize(x_0,lower_level_problem,upper_level_problem,param);
  s.algorithm = @(x_0,lower_level_problem,upper_level_problem,sol,s,param) nonsmooth_trust_region_algorithm(lower_level_problem,upper_level_problem,sol,s,param);
  s.finalize = @(x_0,lower_level_problem,upper_level_problem,sol,s,param) sol;
end

function [sol,s,param] = nonsmooth_trust_region_initialize(x_0,lower_level_problem,upper_level_problem,param)

  s.x_n = {};
  sol = x_0;
  s.radius = param.radius;
  s.hess = 0;
  s.yold = lower_level_problem.solve(sol);
  s.gradold = upper_level_problem.adjoint(s.yold,sol);

  % Test if the min radius is defined
  if ~isfield(param,'minradius')
    param.minradius = 1e-4;
  end

  % Test if lower level problem has a solve method
  if ~isfield(lower_level_problem, 'solve')
    error('Lower Level Problem struct does not provide a SOLVE method.')
  end

  % Test if upper level problem has an adjoint method
  if ~isfield(upper_level_problem, 'adjoint')
    error('Upper Level Problem struct does not provide an ADJOINT method.')
  end

  % Test if upper level problem has a slack method
  if ~isfield(upper_level_problem, 'slack')
    error('Upper Level Problem struct does not provide an SLACK method.')
  end
end

function [sol,s] = nonsmooth_trust_region_algorithm(lower_level_problem,upper_level_problem,sol,s,param)

  % Solving the state equation (lower level solver)
  y = lower_level_problem.solve(sol);

  % Solving the adjoint state
  s.grad = upper_level_problem.adjoint(y,sol);

  % Getting current cost
  cost = upper_level_problem.eval(y,sol);

  if s.radius >= param.minradius
    % Get BFGS Matrix
    % dy = y-s.yold;
    % %tt = hess*dy;
    % dgrad = s.grad-s.gradold;
    % if norm(y-s.yold) ~= 0
    %   s.hess = bfgs(s.hess,dgrad,dy);
    % end

    % Trust Region Step Calculation (Solving TR Subproblem)
    step = tr_subproblem(s.grad,s.hess,s.radius);

    % Record previous step
    s.yold = y;
    s.gradold = s.grad;

    % Trust Region Modification
    pred = -s.grad'*step-0.5*step'*s.hess*step;
    next_y = lower_level_problem.solve(sol+step);
    next_cost = upper_level_problem.eval(next_y,sol+step);
    ared = cost-next_cost;
    rho = ared/pred;

    % Change size of the region
    if rho > param.eta2
      sol = sol + step;
      s.radius = param.gamma2*s.radius;
    elseif rho <= param.eta1
      s.radius = param.gamma1*s.radius;
    else
      %sol = sol + step;
      s.radius = param.gamma1*s.radius;
    end

    fprintf('sol = %f, grad = %f, radius = %f, rho = %f, step = %f, y = %f, q = %f\n',sol,s.grad,s.radius,rho,step,y,upper_level_problem.slack(y,sol));

  else
    grad = [0.4*sol;0.5*(y-1)+0.4*sol];

    [xi,step] = tr_subproblem_complex(grad,s.radius);
    stationarity = tr_complex_stationarity_measure(grad);

    % Record previous step
    s.yold = y;
    s.gradold = s.grad;

    % Trust Region Modification
    pred = -xi;
    next_y = lower_level_problem.solve(sol+step);
    next_cost = upper_level_problem.eval(next_y,sol+step);
    ared = cost-next_cost;
    rho = ared/pred;

    % Change size of the region
    if rho > param.eta2
      sol = sol + step;
      s.radius = param.gamma2*s.radius;
    elseif rho <= param.eta1
      s.radius = param.gamma1*s.radius;
    else
      %sol = sol + step;
      s.radius = param.gamma1*s.radius;
    end

    fprintf('sol = %f, grad = %f, radius = %f, rho = %f, step = %f, stat = %f, y = %f , q = %f\n',sol,norm(grad),s.radius,rho,step,stationarity,y,upper_level_problem.slack(y,sol));
  end

end

function step = tr_subproblem(grad,hess,radius)
  % Step calculation
  sn = -hess\grad;
  predn = -grad'*sn-0.5*sn'*hess*sn;
  if grad'*hess*grad <= 0
    t = radius/norm(grad);
  else
    t = min(norm(grad).^2/(grad'*hess*grad),radius/(norm(grad)));
  end
  sc = -t*grad;
  predc = -grad'*sc-0.5*sc'*hess*sc;

  % Step Selection
  if norm(sn)<=radius && predn >= 0.8*predc
    step = sn;
  else
    step = sc;
  end
end

function B = bfgs(B,y,s)
  alpha = 1/(y'*s);
  beta = 1/(s'*B*s);
  u = y*y';
  v = (B*s)*(s'*B');
  B = B + alpha*u - beta*v;
end

function A = constraints_matrix(grad)
  [m,n] = size(grad);
  A1 = [-ones(m,1) grad];
  A2 = [zeros(n,1) eye(n)];
  A3 = [zeros(n,1) -eye(n)];
  A = [A1;A2;A3];
end

function [xi,step] = tr_subproblem_complex(grad,radius)
  options = optimset('Display','none');
  [m,n] = size(grad);
  c = [1;zeros(n,1)];
  b = [zeros(m,1);radius*ones(2*n,1)];
  A = constraints_matrix(grad);
  linsol = linprog(c,A,b,[],[],[],[],[],options);
  xi = linsol(1);
  step = linsol(2:end);
end

function xi = tr_complex_stationarity_measure(grad)
  [xi, ~] = tr_subproblem_complex(grad,1);
  xi = -xi;
end

% function nXi = xi(p,m,n)
%   p = reshape(p,m*n,2);
%   a = sqrt(sum(p.^2,2));
%   nXi = [a;a];
% end
%
% function prod = outer_product(p,q,m,n)
%   p = reshape(p,m*n,2);
%   q = reshape(q,m*n,2);
%   a = p(:,1).*q(:,1);
%   b = p(:,1).*q(:,2);
%   c = p(:,2).*q(:,1);
%   d = p(:,2).*q(:,2);
%   prod = [spdiags(a,0,m*n,m*n) spdiags(b,0,m*n,m*n); spdiags(c,0,m*n,m*n) spdiags(d,0,m*n,m*n)];
% end
%
% function [adj,grad] = adjoint_solver(u,original,sol)
%   % Get the adjoint state
%   [m,n] = size(u);
%   nabla = gradient_matrix(m,n);
%   Ku = nabla*u(:);
%   nKu = xi(Ku,m,n);
%   act = (nKu<1e-7);
%   inact = 1-act;
%   Act = spdiags(act,0,2*m*n,2*m*n);
%   Inact = spdiags(inact,0,2*m*n,2*m*n);
%   denominador = Inact*nKu+act;
%   prodKuKu = outer_product(Ku./(denominador.^3),Ku,m,n);
%   A = speye(m*n);
%   B = nabla';
%   C = sol*Inact*(prodKuKu-spdiags(1./denominador,0,2*m*n,2*m*n))*nabla;
%   D = speye(2*m*n);
%   E = Act*nabla;
%   F = sparse(2*m*n,2*m*n);
%   Adj = [A B;C D;E F];
%   Track = [u(:)-original(:);sparse(4*m*n,1)];
%   mult = Adj\Track;
%   adj = mult(1:n*m);
%   % Calculating the gradient
%   Kp = nabla*adj;
%   aux = Inact*(Ku./denominador);
%   grad = sol*0.001*m*n - aux'*Kp;
% end
