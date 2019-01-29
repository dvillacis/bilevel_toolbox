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
  s.res = 1;

  % Test if the min radius is defined
  if ~isfield(param,'minradius')
    param.minradius = 1e-4;
  end

  % Test if lower level problem has a solve method
  if ~isfield(lower_level_problem, 'solve')
    error('Lower Level Problem struct does not provide a SOLVE method.')
  end

  % Test if upper level problem has an gradient method
  if ~isfield(upper_level_problem, 'gradient')
    error('Upper Level Problem struct does not provide an GRADIENT method.')
  end

end

function [sol,s] = nonsmooth_trust_region_algorithm(lower_level_problem,upper_level_problem,sol,s,param)

    % Solving the state equation (lower level solver)
    u = lower_level_problem.solve(sol);
    u = u(:);

    % Solving the gradient
    s.grad = upper_level_problem.gradient(u,sol,s.radius);

    % Getting current cost
    cost = upper_level_problem.eval(u,sol);

    if s.radius >= param.minradius

        % Hessian Matrix approximation
        if isfield(s,'solprev') && norm(s.hess)~= 0
            dsol = sol - s.solprev;
            r = s.grad-s.gradprev;
            if s.solprev - sol ~= 0 && dsol'*r > 0
                s.hess = bfgs(s.hess,dsol,r);
            end
        end

        % Trust Region Step Calculation (Solving TR Subproblem)
        step = tr_subproblem(s.grad,s.hess,s.radius);

        % Project the step to the positive quadrant
        step(sol+step < 0) = 0;

        % Record previous step
        s.solprev = sol;
        s.gradprev = s.grad;

        % Trust Region Modification
        pred = -s.grad'*step-0.5*step'*s.hess*step;
        next_u = lower_level_problem.solve(sol+step);
        next_cost = upper_level_problem.eval(next_u,sol+step);
        ared = cost-next_cost;
        rho = ared/pred;
        
        if size(sol,1)>1 || size(sol,2)>1
            fprintf('l2_cost = %f, norm_sol = %f, norm_grad = %f, radius = %f, rho = %f\n',cost,norm(sol),norm(s.grad),s.radius,rho);
        else
            fprintf('l2_cost = %f, sol = %f, grad = %f, radius = %f, rho = %f\n',cost,sol,s.grad,s.radius,rho);
        end

        % Change size of the region
        if rho > param.eta2
          sol = sol + step;
          s.radius = param.gamma2*s.radius;
          s.current_l2_cost = next_cost;
        elseif rho <= param.eta1
          s.radius = param.gamma1*s.radius;
          s.current_l2_cost = cost;
        else
          %sol = sol + step;
          s.radius = param.gamma1*s.radius;
          s.current_l2_cost = cost;
        end

    else

        [xi,step] = tr_subproblem_complex(s.grad,s.hess,s.radius);
        psi = tr_complex_stationarity_measure(s.grad,s.hess);

        % Trust Region Modification
        pred = -xi-0.5*step'*s.hess*step;
        next_y = lower_level_problem.solve(sol+step);
        next_cost = upper_level_problem.eval(next_y,sol+step);
        ared = cost-next_cost;

        if (psi > norm(s.grad)*s.radius)
            rho = ared/pred;
        else
            rho = 0;
        end

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

        fprintf('COMPLEX: sol = %f, grad = %f, radius = %f, rho = %f, step = %f, psi=%f\n',sol,norm(s.grad),s.radius,rho,norm(step),psi);
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

function [c,ceq] = norm_constraint(x,radius)
    c = norm(x)-radius;
    ceq = [];
end

function [xi,step] = tr_subproblem_complex(grad,hess,radius)
    [m,n] = size(grad);
    B = spdiags(ones(n+1,1),1,m,n+1);
    H = B'*hess*B;
    f = [1;zeros(n,1)];
    b = zeros(m,1);
    A = [-ones(m,1),grad];
    obj = @(x) 0.5*x'*H*x + f'*x;
    nonloc = @(x) norm_constraint(x,radius);
    options = optimoptions('fmincon','Display','none');
    [x,fval] = fmincon(obj,[0;0],A,b,[],[],[],[],nonloc,options);
    xi = x(1);
    step = x(2:end);
end

function xi = tr_complex_stationarity_measure(grad,hess)
  [xi, ~] = tr_subproblem_complex(grad,0,1);
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
