function [G,Gx,Gy] = gradient_matrix(M,N)

  %Set up N-D finite (forward) difference matrix with Neumann boundary conditions
  %% image
  np = M*N;
  idx = reshape(1:np, M,N);

  %% u_{i,j+1} - u_{i,j}
  idx1 = idx;
  idx2 = idx(:,[2:N,N]);

  Gx = -sparse(1:np, idx1(:), ones(np,1), np, np) + ...
        sparse(1:np, idx2(:), ones(np,1), np, np);

  %% u_{i+1,j} - u_{i,j}
  idx1 = idx;
  idx2 = idx([2:M,M],:);

  Gy = -sparse(1:np, idx1(:), ones(np,1), np, np) + ...
        sparse(1:np, idx2(:), ones(np,1), np, np);

  G = [Gx;Gy];

end
