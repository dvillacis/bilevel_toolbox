function res=dynsel(fast, normal, varargin)
    if isa(varargin{1}, 'cvx')
        res=normal(varargin{:});
    else
        res=fast(varargin{:});
    end
end
