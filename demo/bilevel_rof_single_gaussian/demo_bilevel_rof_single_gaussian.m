clear all;
close all;
clc;

init_bilevel_toolbox();

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset_single_gaussian','*_circle_original.png','*_circle_noisy.png');

%% Load input image
original = dataset.get_target(1);
noisy = dataset.get_corrupt(1);

% Define lower level problem
lower_level_problem.solve = @(alpha) solve_lower_level(alpha,noisy);

% Define upper level problem
upper_level_problem.eval = @(u,alpha) 0.5*norm(u(:)-original(:)).^2;
upper_level_problem.gradient = @(u,alpha,radius) solve_gradient(u,alpha,radius,original);

%% Solving the bilevel problem
bilevel_param.verbose = 2;
bilevel_param.maxit = 1000;
bilevel_param.tol = 1e-4;
bilevel_param.algo = 'NONSMOOTH_TRUST_REGION';
bilevel_param.radius = 1;
bilevel_param.minradius = 0.0001;
bilevel_param.gamma1 = 0.5;
bilevel_param.gamma2 = 1.5;
bilevel_param.eta1 = 0.01;
bilevel_param.eta2 = 0.93;
alpha = 0.1;
[sol,info] = solve_bilevel(alpha,lower_level_problem,upper_level_problem,bilevel_param);

optimal_sol = solve_lower_level(sol,noisy);
%% Plotting the solution
figure(1)
subplot(1,3,1)
imagesc_gray(original,1,'Original Image');
subplot(1,3,2)
imagesc_gray(noisy,1,'Noisy Image');
subplot(1,3,3)
imagesc_gray(optimal_sol,1,'Denoised Image');

%% Auxiliary functions

function y = solve_lower_level(alpha,noisy)
  param_lower_level.maxit = 2000;
  param_lower_level.alpha = alpha;
  param_lower_level.verbose = 0;
  y = solve_rof_cp_single_gaussian(noisy,param_lower_level);
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

function grad = solve_gradient(u,alpha,~,original)
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
  C = alpha*Inact*(prodKuKu-spdiags(1./denominador,0,2*m*n,2*m*n))*nabla;
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
  grad = alpha*0.001*m*n - aux'*Kp;
end
