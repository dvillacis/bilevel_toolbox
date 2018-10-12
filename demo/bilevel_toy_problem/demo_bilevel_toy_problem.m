clear all;
close all;
clc;

init_bilevel_toolbox();

% Define lower level problem
lower_level_problem.solve = @(u) solve_lower_level(u);

% Define upper level problem
zd = 1;
alpha = 0.4;
A = 2;
upper_level_problem.eval = @(y,u) 0.5*norm(y-zd).^2 + 0.5*alpha*norm(u).^2;
upper_level_problem.adjoint = @(y,u,radius) solve_adjoint(y,u,radius,zd,alpha,A);
upper_level_problem.slack = @(y,u) u-A*y;

% Initial control
u = -5;

% Plotting the bilevel cost function
lb = -3;
ub = 3;
c = [];
for uu = lb:0.01:ub
  yy = lower_level_problem.solve(uu);
  cc = upper_level_problem.eval(yy,uu);
  c = [c cc];
end
plot(lb:0.01:ub,c);
hold on;

% Define bilevel parameters
bilevel_param.verbose = 2;
bilevel_param.maxit = 1000;
bilevel_param.tol = 1e-4;
bilevel_param.algo = 'NONSMOOTH_TRUST_REGION';
bilevel_param.radius = 1;
bilevel_param.minradius = 0.1;
bilevel_param.gamma1 = 0.5;
bilevel_param.gamma2 = 1.5;
bilevel_param.eta1 = 0.01;
bilevel_param.eta2 = 0.93;

% Solve the bilevel problem
[sol,info] = solve_bilevel(u,lower_level_problem,upper_level_problem,bilevel_param);
y = lower_level_problem.solve(sol);
plot(sol,upper_level_problem.eval(y,sol),'r*');


% Auxiliary functions

% Lower Level Solver
function y = solve_lower_level(u)
  if u >= 1
    y = 0.5*(u-1);
  elseif u <= -1
    y = 0.5*(u+1);
  else
    y = 0;
  end
end

% Adjoint Solver
function grad = solve_adjoint(y,u,radius,zd,alpha,A)
    % Getting the very active and possibly biactive sets
    slack = u-A*y;
    active = find(abs(slack) < 1-0.5*radius);
    biactive = find(abs(y) <= 0.5*radius && abs(slack) >= 1-0.5*radius);
    
    if isempty(active) && isempty(biactive)
        grad = 0.5*(y-zd)+alpha*u;
    elseif ~isempty(active) && isempty(biactive)
        grad = alpha*u;
    else
        grad = [0.5*(y-zd)+alpha*u;alpha*u];
    end
end
