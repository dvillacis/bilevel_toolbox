function [y] = l2_norm(x)
%L2_NORM L2 norm of a vector
%   Detailed explanation goes here

% Check if the size of the vector is even
if ~isvector(x)
  error('x must be a vector.')
end

[N,M] = size(x);
if mod(N,2) ~= 0
  error('x has odd size.')
end

x_ = reshape(x,N/2,2);
y = sqrt(sum(x_.^2,2));

end
