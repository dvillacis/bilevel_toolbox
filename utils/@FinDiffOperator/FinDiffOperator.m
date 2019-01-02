classdef FinDiffOperator < MatrixOperator
    %FINDIFFOPERATOR Finite Differences Operator
    %   Finite differences matrix operator

    properties
        Dim
        Method
        Weight
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

        function op_val = val(obj,x)
            %EVAL Evaluation of the operator
            %   Detailed explanation goes here
            if useopt && method(1)=='f' && method(2)=='n'
                op_val = obj.selmex(usemex, cvxversion, weight*fastdiff(x, 'gfn'), weight*fwddiff(x));
            elseif useopt && method(1)=='b' && method(2)=='n'
                op_val = weight*bwddiff(x);
            elseif useopt && method(1)=='c' && method(2)=='n'
                op_val = weight*ctrdiff(x);
            else
                dualdim=[dim, 2];
                grad=weight*diff2d(dim, method);
                op_val=weight*reshape(grad*x(:), dualdim);
            end

        end

        function op_conj = conj(obj,y)
            %EVAL_CONJ Evaluation of the conjugate of the operator
            %   Detailed explanation goes here
            if useopt && method(1)=='f' && method(2)=='n'
                op_conj=obj.selmex(usemex, cvxversion, weight*fastdiff(y, '*fn'), weight*fwddiff_conj(y));
            elseif useopt && method(1)=='b' && method(2)=='n'
                op_conj = weight*bwddiff_conj(y);
            elseif useopt && method(1)=='c' && method(2)=='n'
                op_conj = weight*ctrdiff_conj(y);
            else
                grad = weight*diff2d(dim, method);
                op_conj = weight*reshape(grad'*y(:), dim);
            end
        end

        function grad = matrix(obj)
            %MATRIX Return a matrix representation of the operator
            %   Detailed explanation goes here
            grad = diff2d(obj.Dim,obj.Method);
        end

        function bnd = bound(obj)
            bnd = norm(obj.matrix());
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
