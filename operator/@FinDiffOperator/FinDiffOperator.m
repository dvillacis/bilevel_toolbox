classdef FinDiffOperator < Operator
    %FINDIFFOPERATOR Finite Differences Operator
    %   Finite differences matrix operator

    properties
        Dim
        Method
        Weight
        UseOpt
        UseMex
        CVXVersion
        %Grad
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
            obj.Weight = 1;
            %obj.Grad = gradient_matrix(dim(1),dim(2));
        end

        function res = val(obj,x)
            %EVAL Evaluation of the operator
            %   Detailed explanation goes here
            if obj.UseOpt && obj.Method(1)=='f' && obj.Method(2)=='n'
                op_val = obj.selmex(obj.UseMex, obj.CVXVersion, @(x)fastdiff(x, 'gfn'), @(x) obj.fwddiff(x));
                res = op_val(x);
            elseif obj.UseOpt && obj.Method(1)=='b' && obj.Method(2)=='n'
                res = obj.Weight*bwddiff(x);
            elseif obj.UseOpt && obj.Method(1)=='c' && obj.Method(2)=='n'
                res = obj.Weight*ctrdiff(x);
            else
                error("Will use unoptimized method- deprecated");
%             elseif obj.UseOpt && obj.Method(1)=='c' && obj.Method(2)=='n'
%                 res = obj.Weight*ctrdiff(x);
%             else
%                 dualdim=[obj.Dim, 2];
%                 %grad=obj.Weight*diff2d(obj.Dim, obj.Method);
%                 res = obj.Weight*reshape(obj.Grad*x(:), dualdim);
%                 %res = obj.Weight*reshape(grad*x(:), dualdim);
            end

        end

        function res = conj(obj,y)
            %EVAL_CONJ Evaluation of the conjugate of the operator
            %   Detailed explanation goes here
            if obj.UseOpt && obj.Method(1)=='f' && obj.Method(2)=='n'
                op_conj=obj.selmex(obj.UseMex, obj.CVXVersion, @(y) fastdiff(y, '*fn'), @(y) obj.fwddiff_conj(y));
                res = op_conj(y);
            elseif obj.UseOpt && obj.Method(1)=='b' && obj.Method(2)=='n'
                res = obj.Weight*bwddiff_conj(y);
            elseif obj.UseOpt && obj.Method(1)=='c' && obj.Method(2)=='n'
                res = obj.Weight*ctrdiff_conj(y);
            else
                error("Will use unoptimized method- deprecated");
%             elseif obj.UseOpt && obj.Method(1)=='c' && obj.Method(2)=='n'
%                 res = obj.Weight*ctrdiff_conj(y);
%             else
%                 %grad = obj.Weight*diff2d(obj.Dim, obj.Method);
%                 res = obj.Weight*reshape(obj.Grad'*y(:), obj.Dim);
%                 %res = obj.Weight*reshape(grad'*y(:), obj.Dim);
%             end
            end
        end

%         function grad = matrix(obj)
%             %MATRIX Return a matrix representation of the operator
%             %   Detailed explanation goes here
%             grad = diff2d(obj.Dim,obj.Method);
%         end

        function bnd = bound(obj)
            %bnd = norm(obj.matrix());
            bnd = sqrt(8);
        end
    end

    methods(Static)

        function fn = selmex(usemex,cvxversion,fast,normal)
            %SELMEX Solver selection
            %   Detailed explanation goes here
            if ~usemex
                fn = normal;
            elseif ~cvxversion
                fn = fast;
            else
                fn=@(x) dynsel(fast,normal,x);
            end

        end

        function res = fwddiff(x)
            res=zeros([size(x), 2]);
            if isa(x, 'cvx')
                res=cvx(res);
            end
            res(1:(end-1), :, 1)=x(2:end, :)-x(1:(end-1), :);
            res(:, 1:(end-1), 2)=x(:, 2:end)-x(:, 1:(end-1));
        end

        function res=fwddiff_conj(y)
            sz=size(y);
            assert(sz(3)==2);
            res=zeros(sz(1), sz(2));
            if isa(y, 'cvx')
                res=cvx(res);
            end
            res(1:(end-1), :)=                  -y(1:(end-1), :, 1);
            res(2:end, :)    =res(2:end, :)     +y(1:(end-1), :, 1);
            res(:, 1:(end-1))=res(:, 1:(end-1)) -y(:, 1:(end-1), 2);
            res(:, 2:end)    =res(:, 2:end)     +y(:, 1:(end-1), 2);
        end

        function res=bwddiff(x)
            res=zeros([size(x), 2]);
            if isa(x, 'cvx')
                res=cvx(res);
            end
            res(2:end, :, 1)=x(2:end, :)-x(1:(end-1), :);
            res(:, 2:end, 2)=x(:, 2:end)-x(:, 1:(end-1));
        end

        function res=bwddiff_conj(y)
            sz=size(y);
            assert(sz(3)==2);
            res=zeros(sz(1), sz(2));
            if isa(y, 'cvx')
                res=cvx(res);
            end

            res(2:end, :)    =                  +y(2:end, :, 1);
            res(1:(end-1), :)=res(1:(end-1), :) -y(2:end, :, 1);
            res(:, 2:end)    =res(:, 2:end)     +y(:, 2:end, 2);
            res(:, 1:(end-1))=res(:, 1:(end-1)) -y(:, 2:end, 2);
        end

        function res=ctrdiff(x)
            res=zeros([size(x), 2]);
            if isa(x, 'cvx')
                res=cvx(res);
            end

            res(2:(end-1), :, 1)=(x(3:end, :)-x(1:(end-2), :))/2;
            res(:, 2:(end-1), 2)=(x(:, 3:end)-x(:, 1:(end-2)))/2;
        end

        function res=ctrdiff_conj(y)
            sz=size(y);
            assert(sz(3)==2);
            res=zeros(sz(1), sz(2));
            if isa(y, 'cvx')
                res=cvx(res);
            end

            res(3:end, :)    =                  +y(2:(end-1), :, 1)/2;
            res(1:(end-2), :)=res(1:(end-2), :) -y(2:(end-1), :, 1)/2;
            res(:, 3:end)    =res(:, 3:end)     +y(:, 2:(end-1), 2)/2;
            res(:, 1:(end-2))=res(:, 1:(end-2)) -y(:, 2:(end-1), 2)/2;
        end
    end
end
