function res=norm2sq(x)
    fn=@(x) dynsel(@stablenormsq, @basicnormsq, x);
    if isreal(x)
        res=fn(x);
    else
        res=fn(real(x))+fn(imag(x));
    end
end

function res=basicnormsq(x)
    res=sum(x(:).*x(:));
end