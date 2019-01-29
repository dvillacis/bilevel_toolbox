clear all;
close all;
clc;

init_bilevel_toolbox();

%% Load dataset
dataset = DatasetInFolder('data/circle_dataset_single_gaussian','*_circle_original.png','*_circle_noisy.png');
%dataset = DatasetInFolder('data/smiley','*_smiley_original.png','*_smiley_noisy.png');

%% Load input image
original = dataset.get_target(1);
noisy = dataset.get_corrupt(1);

% Define lower level problem
lower_level_problem.solve = @(lambda) solve_lower_level(lambda,noisy);

% Define upper level problem
upper_level_problem.eval = @(u,lambda) 0.5*norm(u(:)-original(:)).^2 + 0.0001*norm(lambda);
upper_level_problem.gradient = @(u,lambda,params) solve_gradient(u,lambda,original,noisy,params);
upper_level_problem.dataset = dataset;

%% Solving the bilevel problem
bilevel_param.verbose = 2;
bilevel_param.maxit = 100;
bilevel_param.tol = 1e-4;
bilevel_param.algo = 'NONSMOOTH_TRUST_REGION';
bilevel_param.radius = 0.01;
bilevel_param.minradius = 0.00001;
bilevel_param.gamma1 = 0.5;
bilevel_param.gamma2 = 1.5;
bilevel_param.eta1 = 0.01;
bilevel_param.eta2 = 0.80;
lambda = 0.1;
[sol,info] = solve_bilevel(lambda,lower_level_problem,upper_level_problem,bilevel_param);

optimal_sol = solve_lower_level(sol,noisy);
%% Plotting the solution
[M,N]=size(original);
optimal_sol = reshape(optimal_sol,M,N);
figure(1)
subplot(1,3,1)
imagesc_gray(original,1,'Original Image');
subplot(1,3,2)
imagesc_gray(noisy,1,'Noisy Image');
subplot(1,3,3)
imagesc_gray(optimal_sol,1,'Denoised Image');

%% Auxiliary functions

function y = solve_lower_level(lambda,noisy)
    param_solver.verbose = 0;
    param_solver.maxiter = 2000;
    %% Define the cell matrices
    [M,N] = size(noisy);
    K = speye(M*N);
    z = noisy(:);
    gradient = FinDiffOperator([M,N],'fn');
    B = gradient.matrix();
    q = zeros(2*M*N,1);
    alpha = 1;
    gamma = 0;
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

function adj = calculate_adjoint(u,Ku,nKu,nabla,lambda,original,Inact,Act,act,M,N)

    denominador = Inact*nKu+act;
    prodKuKu = outer_product(Ku./(denominador.^3),Ku,M,N);
    A = lambda*speye(M*N);
    B = nabla';
    C = -Inact*(prodKuKu-spdiags(1./denominador,0,2*M*N,2*M*N))*nabla;
    D = speye(2*M*N);
    E = Act*nabla;
    F = sparse(2*M*N,2*M*N);
    Adj = [A B;C D;E F];
    Track = [u(:)-original(:);sparse(4*M*N,1)];
    mult = Adj\Track;
    adj = mult(1:N*M);
end

function grad = solve_gradient(u,lambda,original,noisy,params)

    [M,N] = size(noisy);
    nabla = gradient_matrix(M,N);
    Ku = nabla*u(:);
    nKu = xi(Ku,M,N);

    if params.complex_model == false
        act = (nKu<1e-7);
        inact = 1-act;
        Act = spdiags(act,0,2*M*N,2*M*N);
        Inact = spdiags(inact,0,2*M*N,2*M*N);

        % Get the adjoint state
        adj = calculate_adjoint(u,Ku,nKu,nabla,lambda,original,Inact,Act,act,M,N);
    else
        act = (nKu<1e-7);
        inact = 1-act;
        Act = spdiags(act,0,2*M*N,2*M*N);
        Inact = spdiags(inact,0,2*M*N,2*M*N);

        % Get the adjoint state
        adj = calculate_adjoint(u,Ku,nKu,nabla,lambda,original,Inact,Act,act,M,N);
    end

    % Calculating the gradient
    beta = 0.0001;
    grad = (noisy(:)-u(:))'*adj + beta*lambda;
end
