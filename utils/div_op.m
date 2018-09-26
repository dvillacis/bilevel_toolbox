function I = div_op(x, dx, dy)
%DIV_OP Divergence operator in 2 dimensions
%   Usage:  I = div_op(dx, dy)
%
%   Input parameters:
%         dx    : Gradient along x
%         dy    : Gradient along y
%
%   Output parameters:
%         I     : Output divergence image
%
%   Compute the 2-dimensional divergence of an image. If a cube is given,
%   it will compute the divergence of all images in the cube.
%
%   Warning: computes the divergence operator defined as minus the adjoint
%   of the gradient
%
%   ..      div  = - grad'
%
%   .. math:: \text{div} = - \nabla^*
%
    [M,N] = size(x);
    G = gradient_matrix(M,N);
    I = G'*[dx;dy];
    I = reshape(I,M,N);

end
