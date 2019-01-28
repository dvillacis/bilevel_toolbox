clear all;
close all;
clc;

init_bilevel_toolbox();

%% Load dataset
%dataset = DatasetInFolder('data/circle_dataset_single_gaussian','*_circle_original.png','*_circle_noisy.png');
dataset = DatasetInFolder('data/smiley','*_smiley_original.png','*_smiley_noisy.png');

%% Load input image
original = dataset.get_target(1);
noisy = dataset.get_corrupt(1);
[M,N]=size(original);

% Define lower level problem
lower_level_problem.solve = @(lambda) solve_lower_level(lambda,noisy);

% Define upper level problem
upper_level_problem.eval = @(u,lambda) 0.5*norm(u(:)-original(:)).^2 + 0.5*0.0001*norm(lambda).^2;
upper_level_problem.gradient = @(u,lambda,radius) solve_gradient(u,lambda,radius,original,noisy);
upper_level_problem.dataset = dataset;

%% Solving the bilevel problem
bilevel_param.verbose = 2;
bilevel_param.maxit = 100;
bilevel_param.tol = 1e-2;
bilevel_param.algo = 'NONSMOOTH_TRUST_REGION';
bilevel_param.radius = 0.5;
bilevel_param.minradius = 0.00001;
bilevel_param.gamma1 = 0.5;
bilevel_param.gamma2 = 2.0;
bilevel_param.eta1 = 0.01;
bilevel_param.eta2 = 0.70;
lambda1 = 1*ones(0.5*M*N,1);
lambda2 = 5*ones(0.5*M*N,1);
lambda = vertcat(lambda1,lambda2); % Initial guess
[sol,info] = solve_bilevel(lambda,lower_level_problem,upper_level_problem,bilevel_param);

optimal_sol = solve_lower_level(sol,noisy);
optimal_sol = reshape(optimal_sol,M,N);
%% Plotting the solution
figure(1)
subplot(1,3,1)
imagesc_gray(original,1,'Original Image');
subplot(1,3,2)
imagesc_gray(noisy,1,'Noisy Image');
subplot(1,3,3)
imagesc_gray(optimal_sol,1,'Denoised Image');

figure(2)
[a,b] = meshgrid(1:M,1:N);
surf(a,b,reshape(sol,M,N));

%% Auxiliary functions

function y = solve_lower_level(lambda,noisy)
  %% Solving the Lower Level Problem
  param_solver.verbose = 0;
  param_solver.maxiter = 2000;
  param_solver.tol = 1e-2;

  %% Define the cell matrices
  [M,N] = size(noisy);
  K = speye(M*N);
  z = noisy(:);
  B = gradient_matrix(M,N);
  q = zeros(2*M*N,1);
  alpha = reshape(ones(M,N),M*N,1);
  gamma = 0; % NO Huber regularization

  %% Call the solver
  [y,~] = solve_generic_l1_l2({lambda},{alpha},{K},{B},z,q,gamma,0*noisy(:),param_solver);
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

function grad = solve_gradient(u,lambda,~,original,noisy)
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
  A = spdiags(lambda,0,m*n,m*n);
  B = nabla';
  C = -Inact*(prodKuKu-spdiags(1./denominador,0,2*m*n,2*m*n))*nabla;
  D = speye(2*m*n);
  E = Act*nabla;
  F = sparse(2*m*n,2*m*n);
  Adj = [A B;C D;E F];
  Track = [u(:)-original(:);sparse(4*m*n,1)];
  mult = Adj\Track;
  adj = mult(1:n*m);
  % Calculating the gradient
  grad = (noisy(:)-u(:)).*adj + 0.0001*lambda;
end
