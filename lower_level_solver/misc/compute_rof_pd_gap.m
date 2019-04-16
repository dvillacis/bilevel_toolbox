function [gap, primal, dual] = compute_rof_pd_gap(val_u, u, conj_p, f, alpha)

  norm_val_u = rssq(val_u,3);
  
  primal = alpha * sum(norm_val_u(:)) + 0.5 * sumsqr(u(:)-f(:));
  dual   = f(:)'*conj_p(:) - 0.5 * sumsqr(conj_p(:));
  
  gap = primal-dual;
end
