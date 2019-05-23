classdef PatchOperator < Operator
    %PATCHOPERATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SizeIn
        SizeOut
        KronMatrix
        M
    end
    
    methods
        function obj = PatchOperator(size_in,size_out)
            %PATCHOPERATOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.SizeIn = size_in;
            obj.SizeOut = size_out;
            %assert(obj.SizeIn < obj.SizeOut);
            obj.M = obj.SizeOut./obj.SizeIn;
            obj.KronMatrix = ones(obj.M);
%             N = obj.SizeIn(1);
%             M = obj.SizeOut(1)/N;
%             obj.KronMatrix = ones(M,M);
            
        end
        
        function res = val(obj,x)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            res = kron(x,obj.KronMatrix);
        end
        
        function res = conj(obj,y)
            %res = y(1:size(obj.KronMatrix,1):end,1:size(obj.KronMatrix,2):end)./obj.KronMatrix(1);
            %fun = @(block_struct) mean2(block_struct.data);
            fun = @(block_struct) sum(sum(block_struct.data));
            res = blockproc(full(y),obj.M,fun);
        end
        
        function x_norm = op_norm(~,x)
            x_norm = norm(x);
        end
    end
end

