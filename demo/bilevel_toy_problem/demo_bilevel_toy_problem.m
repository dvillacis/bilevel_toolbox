clear all;
close all;
clc;

init_bilevel_toolbox();

% Define lower level problem
lower_level_problem.solve = @(u) solve_lower_level(u);

% Define upper level problem
upper_level_problem.eval = @(y,u,zd,alpha) 0.5*norm(y-zd).^2 + 0.5*alpha*norm(u).^2;
upper_level_problem.adjoint = @(y,u,zd,alpha) solve_adjoint_upper_level(y,u,zd,alpha);
upper_level_problem.param.zd = 1;
upper_level_problem.param.alpha = 0.4;

% Initial control
u = -10;

% Plotting the bilevel cost function
lb = -3;
ub = 3;
c = [];
for uu = lb:0.01:ub
  yy = lower_level_problem.solve(uu);
  cc = upper_level_problem.eval(yy,uu,upper_level_problem.param.zd,upper_level_problem.param.alpha);
  c = [c cc];
end
plot(lb:0.01:ub,c);
hold on;

% Define bilevel parameters
bilevel_param.verbose = 2;
bilevel_param.maxit = 1000;
bilevel_param.tol = 1e-3;
bilevel_param.algo = 'NONSMOOTH_TRUST_REGION';
bilevel_param.radius = 100;
bilevel_param.minradius = 0.01;
bilevel_param.gamma1 = 0.5;
bilevel_param.gamma2 = 1.5;
bilevel_param.eta1 = 0.01;
bilevel_param.eta2 = 0.99;

% Solve the bilevel problem
[sol,info] = solve_bilevel(u,lower_level_problem,upper_level_problem,bilevel_param);
y = lower_level_problem.solve(sol);
plot(sol,upper_level_problem.eval(y,sol,upper_level_problem.param.zd,upper_level_problem.param.alpha),'r*');


% Auxiliary functions
function y = solve_lower_level(u)
  if u >= 1
    y = 0.5*(u-1);
  elseif u <= -1
    y = 0.5*(u+1);
  else
    y = 0;
  end
end

function grad = solve_adjoint_upper_level(y,u,zd,alpha)
    if u<=1 && u >=-1
        grad = alpha*u;
    else
        grad = 0.5*(y-zd)+alpha*u;
    end
end
