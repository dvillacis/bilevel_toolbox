function [gap, primal, dual] = compute_rof_pd_gap(nabla, u, p, f, lambda, epsilon, M, N)

  div_p = nabla'*p;
  nabla_u = nabla*u;

  nu = rssq(reshape(nabla_u, M*N, 2),2);
  idx = nu < epsilon;
  huber = nu - 0.5*epsilon;
  huber(idx) = (0.5/epsilon) * nu(idx).^2;

  primal = lambda * sum(huber) + 0.5*sumsqr(u-f);
  dual = -0.5*sumsqr(div_p) + f'*div_p - epsilon/(2*lambda*sumsqr(p));

  gap = primal-dual;
end
