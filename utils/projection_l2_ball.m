function [y] = projection_l2_ball(op,x,radius)
%PROJECTION_L2_BALL Project x to a l2 ball of radius
%   Detailed explanation goes here
    den = max(1,op.op_norm(x)./radius);
    y = x./den;
end
