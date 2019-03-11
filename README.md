# Bilevel Toolbox
This is a MATLAB Toolbox designed to test different bilevel optimization problems of the form $ \min_u J(y,u)$

## Installation
To use the toolbox, please execute the script `init_bilevel_toolbox.m` to place all the required scripts into the MATLAB path.

## Dataset
We will use this bilevel optimization techniques in the context of machine learning. In this scenario we use *training data* to define an upper level cost function. In order to provide a helpful interface to this training data, a Dataset class was build and can be instantiated as follows
```matlab
dataset = DatasetInFolder(<training_data_path>,<target_regex>,<data_regex>)
```
This class stores tha path where the data is located, as well as the regular expression required to identify the target and the data.

## Lower Level Problem
In order to define the lower level problem, we need to create a struct that contains a SOLVE method
```matlab
% Define lower level problem
lower_level_problem.solve = @(u) solve_lower_level(u);
```

## Upper Level Problem
This struct must contain the following methods:
* GRADIENT: It returns the gradient for a given parameter.
* EVAL: It calculates the cost function for the upper level problem, this function takes a mandatory solution for the lower level problem.
* DATASET: It is an instance of the Dataset class which specifies the training set to be used in the parameter learning problem.

```matlab
% Define upper level problem
upper_level_problem.eval = @(y,u,zd,alpha) 0.5*norm(y-zd).^2 + 0.5*alpha*norm(u).^2;
upper_level_problem.gradient = @(y,u,zd,alpha) solve_grad_upper_level(y,u,zd,alpha);
upper_level_problem.dataset = dataset;
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
bilevel_param.eta1 = 0.10;
bilevel_param.eta2 = 0.90;

% Solve the bilevel problem
[sol,info] = solve_bilevel(u,lower_level_problem,upper_level_problem,bilevel_param);
```

## Original papers

  1. J.C. de los Reyes and C.-B. Schönlieb, Image denoising: learning the noise model via nonsmooth PDE-constrained optimization, Inverse Problems in Imaging, 7 (4) (2013).
  DOI: [10.3934/ipi.2013.7.1183](http://dx.doi.org/10.3934/ipi.2013.7.1183).

  2. J. C. de Los Reyes, C.-B. Schönlieb and T. Valkonen, Bilevel parameter learning for higher-order total variation regularisation models, Journal of Mathematical Imaging and Vision 57 (2017), 1–25.
  DOI: [10.1007/s10851-016-0662-8](http://dx.doi.org/10.1007/s10851-016-0662-8).

  3. J. C. de Los Reyes, C.-B. Schönlieb and T. Valkonen, The structure of optimal parameters for image restoration problems, Journal of Mathematical Analysis and Applications 434 (2016), 464–500.
  DOI: [10.1016/j.jmaa.2015.09.023](http://dx.doi.org/10.1016/j.jmaa.2015.09.023).

  4. L. Calatroni, C. Cao, J. C. de Los Reyes, C.-B. Schönlieb and T. Valkonen, Bilevel approaches for learning of variational imaging models (2015). Submitted.
  arXiv: [1505.02120](http://arxiv.org/abs/1505.02120).
