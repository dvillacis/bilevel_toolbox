function [gap, primal, dual] = compute_sd_rof_pd_gap(nabla, u, p, f, alpha, M, N)

  div_p = nabla'*p;
  var_u = nabla*u;

  primal = sum(alpha.*sqrt(var_u(1:M*N).^2 + var_u(M*N+1:end).^2)) + 0.5 * sumsqr(u-f);
  dual   = f'*div_p - 0.5 * sumsqr(div_p);

  gap = primal-dual;
end
