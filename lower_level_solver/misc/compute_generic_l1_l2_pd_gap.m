function [gap] = compute_generic_l1_l2_pd_gap(x,y,Ks,Bs,lambda,alpha,z,q)
%COMPUTE_GENERIC_L1_L2_PD_GAP Calculate the primal-dual gap for the
%generic_l1_l2 solver.
%   Detailed explanation goes here
global FUBAR; %% TODO: Remove and rename FUBAR

gamma = 0.01;

% Calculate primal value
primal_1 = 0;
for k=1:length(Ks)
    t = Ks{k}.val(x)-z;
    primal_1 = primal_1 + sum(lambda{k}(:) .* (norm(t(:)).^2));
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
fprintf(' primal = %f, dual = %f\n',primal,dual);
gap = primal+dual;

end
