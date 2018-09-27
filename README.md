# Bilevel Toolbox
This is a MATLAB Toolbox designed to test different bilevel optimization problems of the form $$ \min_u J(y,u)$$

## Installation
To use the toolbox, please execute the script `init_bilevel_toolbox.m` to place all the required scripts into the MATLAB path.

## Lower Level Problem
In order to define the lower level problem, we need to create a struct that contains the following methods: SOLVE and EVAL
```matlab
% Define lower level problem
lower_level_problem.solve = @(u) solve_lower_level(u);
```

## Upper Level Problem
This struct contains the ADJOINT method as well as a EVAL that uses the output from the lower level problem.

```matlab
% Define upper level problem
upper_level_problem.eval = @(y,u,zd,alpha) 0.5*norm(y-zd).^2 + 0.5*alpha*norm(u).^2;
upper_level_problem.adjoint = @(y,u,zd,alpha) solve_adjoint_upper_level(y,u,zd,alpha);
```

## Bilevel Solver
Once both the upper and lower level problems have been properly defined, we can run the bilevel solver. To call this solver some previous parameter configurations are needed.

```matlab
bilevel_param.verbose = 2;
bilevel_param.maxit = 1000;
bilevel_param.tol = 1e-5;
bilevel_param.algo = 'NONSMOOTH_TRUST_REGION';
bilevel_param.radius = 0.5;
bilevel_param.minradius = 0.01;
bilevel_param.gamma1 = 0.5;
bilevel_param.gamma2 = 1.5;
bilevel_param.eta1 = 0.01;
bilevel_param.eta2 = 0.99;

% Solve the bilevel problem
[sol,info] = solve_bilevel(u,lower_level_problem,upper_level_problem,bilevel_param);
```
