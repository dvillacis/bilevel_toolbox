function [ stop ]=test_gamma(gamma)
%TEST_GAMMA test if gamma is correct
%   Usage:  stop = test_gamma(gamma)
%           test_gamma(gamma)
%
%   Input parameters:
%         gamma : number
%   Output parameters:
%         stop  : boolean
%
%   This function test is gamma is stricly positive
%
%   If gamma is negativ, this function return an error. If gamma is zero
%   this function, set stop to 1.
%
%


if gamma<0
    error('gamma can not be negativ!');
end

if gamma==0
    stop = 1;
else
    stop = 0;
end

end
