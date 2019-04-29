function [sol,gap] = solve_generic_l1_l2(lambda,alpha,Ks,Bs,z,q,gamma,xinit,param)
% SOLVE_GENERIC_L1_L2 Generic solver for several image processing problems
% This solver receives an abstract initial structure to support different
% image processing models and solves those by using a Chabolle-Pock algorithm.
% INPUTS
%   lambda: l2 fidelity weights
%   alpha: l1 fidelity weights
%   Ks: cell array of operators corresponding to l2 terms
%   Bs: cell array of operators corresponding to the l1 terms
%   z: l2 data
%   q: l1 data
%   gamma: Huber regularization parameter for l1 terms
%   param: struct with the algorithm specific parameters
%     .maxiter (integer, default: 5000): maximum iteration count
%     .check (integer, default: 100): verbose printing interval
%     .verbose (integer, default: true): verbose flag
%     .tol (integer, default: 1e-3): pd gap stopping criteria
% OUTPUTS
%   sol: minimizer for the optimization problem
%   gap: primal-dual gap values per iteration
%
  global FUBAR; %% TODO: Remove and rename FUBAR to PSEUDO_GAP_STATE
  FUBAR=[];

  % Start the counter
  t1 = tic;

  % Test maxiter parameter
  param = add_default(param, 'maxiter', 5000);

  % Test check parameter
  param = add_default(param, 'check', 100);

  % Test verbose parameter
  param = add_default(param, 'verbose', true);

  % Test tol parameter
  param = add_default(param, 'tol', 1e-3);

  % Test for cell input Ks
  if ~iscell(Ks)
    error('Ks must be a cell array.');
  else
      for k = 1:length(Ks)
          if ~isa(Ks{k},'Operator')
              error('Ks elements must be a Operator class instance.');
          end
      end
  end

  % Test for cell input Bs
  if ~iscell(Bs)
    error('Bs must be a cell array.');
  else
      for k = 1:length(Bs)
          if ~isa(Bs{k},'Operator')
              error('Bs elements must be a Operator class instance.');
          end
      end
  end

  % Test for cell input alpha
  if ~iscell(alpha)
      error('alpha must be a cell array.');
  end

  % Test for cell input alpha
  if ~iscell(lambda)
      error('lambda must be a cell array.');
  end

  % Test for input vector
  if isvector(xinit)
    error('xinit must ve a matrix, not a vector.');
  end

%   L = sqrt(8+1);
%   tau = 0.05/L;
%   sigma = 0.99/(tau*L^2);

  L = sqrt(8);
  tau = 0.01;
  sigma = 1/tau/L^2;

  sol = xinit;
  sol_=sol;

  % Concatenate l2 and l1 matrices
  Kbb = ConcatenatedOperator(Ks{:},Bs{:});
  y = Kbb.val(sol);

  [gap,primal,dual] = compute_generic_l1_l2_pd_gap(sol,y,Ks,Bs,lambda,alpha,z,q);

  if param.verbose > 1
    %fprintf('generic_l1_l2: iter = %4d, primal = %f, dual = %f, gap = %f\n', 0,primal,dual,gap);
    fprintf('generic_l1_l2: iter = %4d, gap = %f\n', 0,gap);
  end

  finished = 0;

  for k = 1:param.maxiter

    % Dual update
    y = y + sigma*Kbb.val(sol_);
    y = calc_prox(y,z,q,Ks,Bs,lambda,alpha,sigma);

    % Primal update
    sol_ = sol;
    sol = sol - tau * Kbb.conj(y);

    % Interpolation step
    sol_ = 2*sol - sol_;

    [ga,~,~] = compute_generic_l1_l2_pd_gap(sol,y,Ks,Bs,lambda,alpha,z,q);
    gap = [gap, ga];

    if mod(k, param.check) == 0 && param.verbose > 1
      %fprintf('generic_l1_l2: iter = %4d, primal = %f, dual = %f, pseudo-gap = %f\n', k, primal, dual, ga);
        fprintf('generic_l1_l2: iter = %4d, pseudo-gap = %f\n', k, ga);
    end

    % Stopping criteria
    if ga < param.tol
        %fprintf('generic_l1_l2: ga = %f, stopping criteria met.\n',ga);
        % Print summary
        if param.verbose>0
            fprintf(['\n ','GENERIC_L1_L2_CHAMBOLLE_POCK',':\n']);
            fprintf(' %i iterations\n', k);
            fprintf(' Primal-Dual Pseudo-Gap: %f \n', gap(end));
            fprintf(' Execution Time: %f \n\n', toc(t1));
        end
        finished = 1;
        break;
    end

  end

  % Print summary
  if param.verbose>0 && finished == 0
    fprintf(['\n ','GENERIC_L1_L2_CHAMBOLLE_POCK',':\n']);
    fprintf(' %i iterations\n', k);
    fprintf(' Primal-Dual Pseudo-Gap: %f \n', gap(end));
    fprintf(' Execution Time: %f \n\n', toc(t1));
  end

end

function prox = calc_prox(y,z,q,Ks,Bs,lambda,alpha,tau)
    for k=1:length(Ks)
        y.elements{k} = (y.elements{k}-tau.*z)./(1+0.5*tau*(1./lambda{k}));% l2 proximal
    end
    for l = 1:length(Bs)
        y.elements{l+length(Ks)} = projection_l2_ball(y.elements{l+length(Ks)}-tau*q,alpha{l});
    end
    prox = y;
end

function [gap,primal,dual] = compute_generic_l1_l2_pd_gap(x,y,Ks,Bs,lambda,alpha,z,q)
%COMPUTE_GENERIC_L1_L2_PD_GAP Calculate the primal-dual gap for the
%generic_l1_l2 solver.
%   Detailed explanation goes here
global FUBAR; %% TODO: Remove and rename FUBAR

gamma = 0.01;

% Calculate primal value
primal_1 = 0;
for k=1:length(Ks)
    t = Ks{k}.val(x)-z;
    primal_1 = primal_1 + sum(lambda{k}(:) .* (t(:).^2));
end
primal_2 = 0;
for l = 1:length(Bs)
    t = Bs{l}.val(x)-q;
    nt = alpha{l}.*l2_norm(t);
    primal_2 = primal_2 + sum(nt(:));
end
primal_reg = 0;
%primal_reg = 0.5*gamma*norm(x).^2;
primal = primal_1 + primal_2 + primal_reg;

% Calculate dual value
dual_1 = 0;
for k=1:length(Ks)
    yk = y.elements{k};
    item = 0.25*((yk.^2)./lambda{k});
    dual_1 = dual_1 + sum(item(:)) + yk(:)'*z(:);
    %dual_1 = dual_1 + 0.25*(1./lambda{k})*norm(yk(:)).^2 + yk(:)'*z(:);
end
dual_2 = 0;
for l = 1:length(Bs)
  yk = y.elements{k+length(Ks)};
  dual_2 = dual_2 + yk(:)' * q(:);
end
Kbb = ConcatenatedOperator(Ks{:},Bs{:});
%dual_reg = 0.5*(1/gamma)*norm2(Kbb.conj(y)).^2;
M=norm(x(:));
if isempty(FUBAR) || FUBAR<M
    FUBAR=M;
end
x_=Kbb.conj(y);
dual_reg = M*norm(x_(:));
dual = dual_1 + dual_2 + dual_reg;

% Calculate the gap
gap = primal+dual;

end
