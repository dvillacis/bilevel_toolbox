
clear;
close all;

verbose = 2; % verbosity level

init_bilevel_toolbox();

%% Creation of the problem
% argmin ||Ax-b||^2 + tau ||x||_1
tau = 1; % regularization parameter

N = 5000; % size of the signal
K = 100; % sparcity level
R = max(4, ceil(log(N))); % constant
fprintf('The compression ratio is: %g\n',N/(R*K));

% Measurements Matrix
A = randn(R * K, N);

% Create a K sparse signal
x = zeros(N,1);
I = randperm(N);
x(I(1:K)) = randn(K,1);
x = x/norm(x);

% Measurements
y = A*x;

%% Defining proximal operators

% Setting the function f2
f2.grad = @(x) 2*A'*(A*x-y);
f2.eval = @(x) norm(A*x-y)^2;
f2.beta = 2*norm(A)^2; % Lipschitz constant of the function

% Setting function f1
param_l1.verbose = verbose - 1;
f1.prox = @(x,T) prox_l1(x, T*tau, param_l1);
f1.eval = @(x) norm(x,1);


%% Solving the problem

% Setting the parameters for the simulation
param_solver.verbose = verbose;
param_solver.maxit = 300;
param_solver.tol = 1e-4;
param_solver.method = 'CHAMBOLLE_POCK';

% Solving the problem
sol = solvep(zeros(N,1), {f1,f2}, param_solver);

% Displaying the result
figure
plot(1:N,x,'o', 1:N,sol,'xr');
legend('Original Signal', 'Reconstructed Signal');
