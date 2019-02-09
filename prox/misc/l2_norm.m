function [y] = l2_norm(x)
%L2_NORM L2 norm of a tensor, componentwise
%   Detailed explanation goes here

x_ = x(:,:,1).^2 + x(:,:,2).^2;
y = sqrt(x_);

end
