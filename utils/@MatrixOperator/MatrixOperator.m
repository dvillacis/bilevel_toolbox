classdef MatrixOperator < Operator
    %MATRIXOPERATOR Operator than can be described with a matrix
    %   Detailed explanation goes here

    properties
        Property1
    end

    methods (Abstract)
        matrix();
    end
end
