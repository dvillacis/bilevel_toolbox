classdef FinDiffOperator < MatrixOperator
    %FINDIFFOPERATOR Finite Differences Operator
    %   Finite differences matrix operator using forward differences and
    %   Newmann boundary conditions.
    
    properties
        Dim
        Method
        UseOpt
        UseMex
        CVXVersion
    end
    
    methods
        function obj = FinDiffOperator(dim,method)
            %FINDIFFOPERATOR Finite Differences Operator Constructor
            %   Detailed explanation goes here
            obj.Dim = dim;
            obj.Method = method;
            obj.UseOpt = true;
            obj.UseMex = true;
            obj.CVXVersion = true;
        end
        
        function outputArg = eval(obj,inputArg)
            %EVAL Evaluation of the operator
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        
        function outputArg = eval_conj(obj,inputArg)
            %EVAL_CONJ Evaluation of the conjugate of the operator
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        
        function outputArg = matrix(obj,inputArg)
            %MATRIX Return a matrix representation of the operator
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

