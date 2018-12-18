classdef MatrixOperator < Operator
    %MATRIXOPERATOR Operator than can be described with a matrix
    %   Detailed explanation goes here

    properties
        
    end

    methods (Abstract)
        matrix();
    end
end
