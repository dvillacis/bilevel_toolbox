function [gap] = compute_generic_l1_l2_pd_gap(x,y,Ks,Bs,lambda,alpha,z,q)
%COMPUTE_GENERIC_L1_L2_PD_GAP Calculate the primal-dual gap for the
%generic_l1_l2 solver.
%   Detailed explanation goes here
global FUBAR; %% TODO: Remove and rename FUBAR

gamma = 0.01;
Kbb = cat(1,Ks{:},Bs{:});

% Calculate primal value
primal_1 = 0;
for k=1:length(Ks)
  primal_1 = primal_1 + lambda .* norm(Ks{k}*x-z).^2;
end
primal_2 = 0;
for l = 1:length(Bs)
  primal_2 = primal_2 + sum(alpha{l} .* l2_norm(Bs{l}*x-q));
end
primal_reg = 0;%0.5*gamma*norm(x).^2;
primal = primal_1 + primal_2 + primal_reg;

% Calculate dual value
index = 0;
dual_1 = 0;
for k=1:length(Ks)
  n = size(Ks{k},1);
  dual_1 = dual_1 + 0.25*(1/lambda)*norm(y(index+1:index+n)).^2 + y(index+1:index+n)' * z;
  index = index+n;
end
dual_2 = 0;
for l = 1:length(Bs)
  n = size(Bs{l},1);
  dual_2 = dual_2 + y(index+1:index+n)' * q;
  index = index+n;
end
%dual_reg = 0.5*(1/gamma)*norm(Kbb'*y).^2;
M=norm(x);
if isempty(FUBAR) || FUBAR<M
    FUBAR=M;
end
dual_reg = M*norm(Kbb'*y);
dual = dual_1 + dual_2 + dual_reg;

% Calculate the gap
gap = primal+dual;
end
