classdef Dataset < handle
    properties
        NumPairs
    end
    methods (Abstract)
        get_pair(id)
        get_target(id)
        get_corrupt(id)
    end
end
