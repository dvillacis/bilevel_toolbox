function [gap, primal, dual] = compute_rof_gap(nabla, u, p, f, lambda, epsilon, M, N)

  div_p = nabla'*p;
  nabla_u = nabla*u;

  nu = sqrt(sum(reshape(nabla_u, M*N, 2).^2,2));
  idx = nu < epsilon;
  huber = nu-epsilon/2;
  huber(idx) = nu(idx).^2/2/epsilon;

  primal = lambda * sum(huber) + 1/2*norm(u-f)^2;
  dual = -norm(div_p)^2/2 + f'*div_p - epsilon/lambda/2*norm(p)^2;

  gap = primal-dual;
end
