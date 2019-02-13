clear all;
close all;
clc;

init_bilevel_toolbox();

%% Load dataset
%dataset = DatasetInFolder('data/circle_dataset_single_gaussian','*_circle_original.png','*_circle_noisy.png');
dataset = DatasetInFolder('data/smiley','*_smiley_original.png','*_smiley_noisy.png');

%% Load input image
% original = dataset.get_target(1);
% noisy = dataset.get_corrupt(1);
[M,N]=size(dataset.get_target(1));

% Define lower level problem
lower_level_problem.solve = @(lambda,dataset) solve_lower_level(lambda,dataset.get_corrupt(1));

% Define upper level problem
upper_level_problem.eval = @(u,lambda,dataset) eval_upper_level(u,lambda,dataset);
upper_level_problem.gradient = @(u,lambda,dataset,params) solve_gradient(u,lambda,dataset,params);
upper_level_problem.dataset = dataset;

%% Solving the bilevel problem
bilevel_param.verbose = 2;
bilevel_param.maxit = 50;
bilevel_param.tol = 1e-2;
bilevel_param.algo = 'NONSMOOTH_TRUST_REGION';
bilevel_param.radius = 5.0;
bilevel_param.minradius = 1.0;
bilevel_param.gamma1 = 0.5;
bilevel_param.gamma2 = 1.5;
bilevel_param.eta1 = 0.10;
bilevel_param.eta2 = 0.60;
lambda = 8.0*triu(ones(M,N))+2.0*tril(ones(M,N)); % Initial guess
%lambda = rand(M,N);
[sol,info] = solve_bilevel(lambda,lower_level_problem,upper_level_problem,bilevel_param);

optimal_sol = solve_lower_level(sol,dataset.get_corrupt(1));

%% Plotting the solution
figure(1)
subplot(1,3,1)
imagesc_gray(dataset.get_target(1),1,'Original Image');
subplot(1,3,2)
imagesc_gray(dataset.get_corrupt(1),1,'Noisy Image');
subplot(1,3,3)
imagesc_gray(optimal_sol,1,'Denoised Image');

figure(2)
[a,b] = meshgrid(1:M,1:N);
surf(a,b,reshape(sol,M,N));

%% Auxiliary functions
function cost = eval_upper_level(u,lambda,dataset)
    original = dataset.get_target(1);
    %cost = 0.5*norm(u(:)-original(:)).^2 + 0.0001*norm(lambda);
    cost = 0.5*norm(u(:)-original(:)).^2 + 0.5*0.1*norm(lambda(:)).^2;
end


function y = solve_lower_level(lambda,noisy)
  %% Solving the Lower Level Problem
  param_solver.verbose = 0;
  param_solver.maxiter = 2000;
  param_solver.tol = 1e-2;

  %% Define the cell matrices
  [M,N] = size(noisy);
  id_op = IdentityOperator([M,N]);
  z = noisy;
  gradient_op = FinDiffOperator([M,N],'fn');
  q = zeros(M,N,2);
  alpha = ones(M,N);
  gamma = 0; % NO Huber regularization

  %% Call the solver
  [y,~] = solve_generic_l1_l2({lambda},{alpha},{id_op},{gradient_op},z,q,gamma,0*noisy,param_solver);
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

function grad = solve_gradient(u,lambda,dataset,params)

    original = dataset.get_target(1);
    noisy = dataset.get_corrupt(1);

    [M,N] = size(noisy);
    nabla = gradient_matrix(M,N);
    A = spdiags(lambda(:),0,M*N,M*N); % Build diagonal matrix with parameters to estimate
    B = nabla';
    Ku = nabla*u(:); %Discrete gradient matrix
    nKu = xi(Ku,M,N); %Discrete euclidean norm

    if params.complex_model == false
        % Get partition active-inactive
        act = (nKu<1e-7); %TODO: Specify a partition of the possible biactive set
        inact = 1-act;
        Act = spdiags(act,0,2*M*N,2*M*N);
        Inact = spdiags(inact,0,2*M*N,2*M*N);

        % Get the adjoint state
        denominador = Inact*nKu+act;
        prodKuKu = outer_product(Ku./(denominador.^3),Ku,M,N);
        C = -Inact*(prodKuKu-spdiags(1./denominador,0,2*M*N,2*M*N))*nabla;
        D = speye(2*M*N);
        E = Act*nabla;
        F = sparse(2*M*N,2*M*N);
        Adj = [A B;C D;E F];
        Track = [u(:)-original(:);sparse(4*M*N,1)];
        mult = Adj\Track;
        adj = mult(1:N*M);
    else
        % Get Active, Strongly Active and Inactive - gamma sets
        gamma = 100;
        act1 = gamma*nKu-1;
        act=spones(max(0,act1(1:M*N)-1/(2*gamma)));
        Act=spdiags(act,0,M*N,M*N);
        inact=spones(min(0,act1(1:M*N)+1/(2*gamma)));
        sact=sparse(1-act-inact);
        Sact=spdiags(sact,0,M*N,M*N);

        % Diagonal matrix corresponding to regularization function
        den=(Act+Sact)*nKu(1:M*N)+inact;
        mk=(act+Sact*(1-gamma/2*(1-gamma*nKu(1:M*N)+1/(2*gamma)).^2))./den;
        Dmi=spdiags(kron(ones(2,1),mk+gamma*inact),0,2*M*N,2*M*N);

        % Negative term in the derivative
        subst=spdiags(act+Sact*(1-gamma/2*(1-gamma*nKu(1:M*N)+1/(2*gamma)).^2),0,M*N,M*N);
        subst=kron(speye(2),subst);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Construction of the Hessian components corresponding to each
        % equation in the optimality system
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Matrix of product (Ku)*(Ku)^T for the Newton step divided by
        % the Frobenius norm to the cube

        H4=outer_product(Ku./kron(ones(2,1),den.^2),Ku,M,N);

        % Matrix of product (Gy)*(Gy)^T for the Newton step

        prodKuKu=outer_product(Ku,Ku,M,N);

        sk2=(Sact*gamma^2*(1-gamma*nKu(1:M*N)+1/(2*gamma)))./(den.^2);
        sk2=spdiags(kron(ones(2,1),sk2),0,2*M*N,2*M*N);

        % Hessian matrix components corresponding to the first equation
        % in the optimality system
        % TODO: This is the second derivative of the l2 norm (check it out!) -> Check derivation!

        hess22=Dmi*nabla-kron(speye(2),(Act+Sact))*Dmi*H4*nabla+sk2*prodKuKu*nabla;

        % Adjoint state is solution of the linear system

        adj=(A+B*hess22)\(u(:)-original(:));

    end

    % Calculating the gradient
    beta = 0.1;
    grad = (noisy-u).*reshape(adj,M,N) + beta*lambda;
end
