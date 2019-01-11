function [H] = bfgs(B,y,s)
%BFGS Summary of this function goes here
%   Detailed explanation goes here
alpha = 1/(y'*s);
beta = 1/(s'*B*s);
u = y*y';
v = (B*s)*(s'*B');
H = B + alpha*u - beta*v;
end

