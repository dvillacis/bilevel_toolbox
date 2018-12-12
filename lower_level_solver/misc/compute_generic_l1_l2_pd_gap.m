function [gap] = compute_generic_l1_l2_pd_gap(x,y,Ks,Bs,lambda,alpha,z,q)
%COMPUTE_GENERIC_L1_L2_PD_GAP Calculate the primal-dual gap for the
%generic_l1_l2 solver.
%   Detailed explanation goes here

% Calculate primal value
primal_1 = sum(lambda.*norm(Ks*x-z).^2);
primal_2 = sum(alpha.*norm(Bs*x-q,1));
primal = primal_1 + primal_2;

% Calculate dual value

gap = primal;
end
