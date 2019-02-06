classdef IdentityOperator < MatrixOperator
    %IDENTITYOPERATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Dim
    end
    
    methods
        function obj = IdentityOperator(dim)
            %IDENTITYOPERATOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.Dim = dim;
        end
        
        function op_val = val(~,x)
            op_val = x;
        end
        
        function op_conj = conj(~,y)
            op_conj = y;
        end
        
        function id = matrix(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            id = speye(obj.Dim(1)*obj.Dim(2));
        end
    end
end

