classdef Operator
    %OPERATOR Abstract class for operators
    %   This abstract class provides an interface for declaring and
    %   handling operators numerically.

    properties

    end

    methods (Abstract)
        val(x);
        conj(x);
    end
end
