classdef Gradient < FinDiffOperator
%GRADIENT Gradient operator based on the finite differences operator

    properties
        ExtraProperty
    end

    methods
        function obj = Gradient(dim,method)
            %GRADIENT Gradient operator constructor
            obj.Dim = dim;
            obj.Method = method;
            obj.Bound = sqrt(8);
            obj.LowerBound = 0;
        end

        function bnd = bound(obj)
            bnd = sqrt(8);
        end
    end
end
