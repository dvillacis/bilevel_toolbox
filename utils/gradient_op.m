function [Gx, Gy] = gradient_op(I)
  %Get the gradients of a given image
  [M,N] = size(I);

  [~,dx,dy] = gradient_matrix(M,N);

  Gx = dx*I(:);
  Gy = dy*I(:);

end
