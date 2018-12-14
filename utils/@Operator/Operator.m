classdef Operator
    %OPERATOR Abstract class for operators
    %   This abstract class provides an interface for declaring and
    %   handling operators numerically.

    properties
        Property1
    end

    methods (Abstract)
        eval(x);
        eval_conj(x);
        bound();
    end
end
