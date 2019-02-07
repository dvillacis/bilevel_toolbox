classdef Tuple

    properties(GetAccess = public, SetAccess = protected)
        elements
        NumElements
    end

    methods(Access = public)
        function obj = Tuple(varargin)
            obj.NumElements = 0;
            if nargin == 0
                error('Tuple Constructor must include at least one element')
            else
                obj.elements = cell([nargin 1]);
                for k = 1:nargin %TODO: Check the element is a matrix or tensor
                    obj.elements{k} = varargin{k};
                    obj.NumElements = obj.NumElements + 1;
                end
            end
        end

        function o3 = plus(o1, o2)
            assert(isa(o1, 'Tuple') && isa(o2, 'Tuple'));
            assert(length(o1)==length(o2));
            sumcell = cell(length(o1.elements),1);
            for k=1:o1.NumElements
                if(size(o1.elements{k})==size(o2.elements{k}))
                    sumcell{k} = o1.elements{k}+o2.elements{k};
                else
                    error('Tuple elements not compatible for adding');
                end
            end
            o3=Tuple(sumcell{:});
        end

        function o3 = minus(o1, o2)
            assert(isa(o1, 'Tuple') && isa(o2, 'Tuple'));
            o3=Tuple(o1.u-o2.u, o1.v-o2.v);
        end

        function o2 = uminus(o1)
            o2=Tuple(-o1.u, -o1.v);
        end

        function o = mtimes(a1, a2)
            if isnumeric(a2)
                o = mtimes(a2, a1);
            else
                assert(isnumeric(a1) && isscalar(a1) && isa(a2, 'Tuple'));
                mtimescell = cell(a2.NumElements,1);
                for k=1:a2.NumElements
                    mtimescell{k} = a1*a2.elements{k};
                end
                o=Tuple(mtimescell{:});
            end
        end

        function o = times(a1, a2)
            if isnumeric(a2)
                o=mtimes(a2, a1);
            elseif isnumeric(a1)
                o=mtimes(a1, a2);
            else
                assert(isa(a1, 'Tuple') && isa(a2, 'Tuple'));
                o=Tuple(a1.u.*a2.u, a1.v.*a2.v);
            end
        end

        function o = mrdivide(a1, a2)
            assert(isnumeric(a2) && isscalar(a2));
            o=Tuple(a1.u/a2, a1.v/a2);
        end

        function n = norm(o)
            n = sqrt(norm(o.u(:))^2+norm(o.v(:))^2);
        end

        function n = norm2(o)
            n = norm(o);
        end

        function n = iprod_real(o1, o2)
            n = iprod_real(o1.u, o2.u)+iprod_real(o1.v, o2.v);
        end

    end
end
