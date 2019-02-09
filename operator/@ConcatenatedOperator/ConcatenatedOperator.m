classdef ConcatenatedOperator
    %CONCATENATEDOPERATOR Summary of this class goes here
    %   Detailed explanation goes here

    properties
        operators
        NumOperators
    end

    methods
        function obj = ConcatenatedOperator(varargin)
            %CONCATENATEDOPERATOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.NumOperators = 0;
            if nargin == 0
                error('ConcatenatedOperator Constructor must include at least one operator');
            else
                obj.operators = cell([nargin 1]);
                for k = 1:nargin
                    if ~isa(varargin{k},'Operator')
                        error('Input argument %d is not an Operator class instance',k);
                    else
                        obj.operators{k} = varargin{k};
                        obj.NumOperators = obj.NumOperators + 1;
                    end
                end
            end
        end

        function op_val = val(obj,x)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            val = cell(obj.NumOperators,1);
            for k = 1:obj.NumOperators
                val{k} = obj.operators{k}.val(x);
            end
            op_val = Tuple(val{:});
        end

        function op_conj = conj(obj,y)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            assert(isa(y,'Tuple'));
            conj = cell(obj.NumOperators,1);
            for k = 1:obj.NumOperators
                conj{k} = obj.operators{k}.conj(y.elements{k});
            end
            conj_tuple = Tuple(conj{:});
            op_conj = conj_tuple.sum_all_elements();
        end

        function ops = get_number_of_operators(obj)
            ops = obj.NumOperators;
        end

    end
end
