function [lambda_out] = patch_operator(lambda_in,size_out)
%PATCH_OPERATOR Extension operator for block constant parameters
%   Detailed explanation goes here
    
    n = size(lambda_in);
    m = size_out./n;

    %assert(n<=size_out);
    
    lambda_out = kron(lambda_in,ones(m));
    


end

function out = inverse_kron(K,array,input_id)

switch input_id
    case 1
        out = K(1:size(K,1)/size(array,1),1:size(K,2)/size(array,2))./array(1);
    case 2
        out = K(1:size(array,1):end,1:size(array,2):end)./array(1);
    otherwise
        error('The Input ID must be either 1 or 2')
end
end

