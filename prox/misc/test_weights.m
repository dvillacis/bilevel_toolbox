function [weights]=test_weights(weights,warning)
%TEST_GAMMA test if the weights are corrects
    if nargin<2
       warning=1;
    end

    if sum(weights<0)
        error('gamma can not be negativ!');
    elseif ~logical(sum(weights(:)~=0)) && warning
        weights=weights+eps;
        fprintf(' WARNING!!! weights is 0. We add eps to weights to keep going...\n');
    end

end
