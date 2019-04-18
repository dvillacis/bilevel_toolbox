function [lambda_out] = patch_operator(lambda_in,size_out)
%PATCH_OPERATOR Extension operator for block constant parameters
%   Detailed explanation goes here

    n = size(lambda_in,1);
    m = size_out/n;
    
    assert(n<=size_out);
    
    lambda_out = kron(lambda_in,ones(m,m));
    


end

