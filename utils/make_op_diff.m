function op=make_op_diff(dim, method, weight)
    assert(length(dim)==2);
    useopt=true;
    usemex=true;
    cvxversion=true;

    if useopt && method(1)=='f' && method(2)=='n'
        op.val=selmex(usemex, cvxversion, @(x) weight*fastdiff(x, 'gfn'), @(x) weight*fwddiff(x));
        op.conj=selmex(usemex, cvxversion, @(x) weight*fastdiff(x, '*fn'), @(x) weight*fwddiff_conj(x));
        op.a_val=selmex(usemex, cvxversion, @(x) weight^2*fastdiff(x, 'afn'), @(x) weight*fwddiff_a(x));
    elseif useopt && method(1)=='b' && method(2)=='n'
        op.val=@(x) weight*bwddiff(x);
        op.conj=@(x) weight*bwddiff_conj(x);
        op.a_val=@(x) weight^2*bwddiff_a(x);
    elseif useopt && method(1)=='c' && method(2)=='n'
        op.val=@(x) weight*ctrdiff(x);
        op.conj=@(x) weight*ctrdiff_conj(x);
        op.a_val=@(x) weight^2*ctrdiff_a(x);
    else
        dualdim=[dim, 2];
        grad=weight*diff2d(dim, method);
        op.val=@(x) weight*reshape(grad*x(:), dualdim);
        op.conj=@(y) weight*reshape(grad'*y(:), dim);
        op.a_val=@(x) (weight^2*4)*reshape(abs(grad)*x(:), dualdim);
    end

    op.bound=weight*sqrt(8);
    op.lowerbound=0;
end

function fn=selmex(usemex, cvxversion, fast, normal)
    if ~usemex
        fn=normal;
    elseif ~cvxversion
        fn=fast;
    else
        fn=@(x) dynsel(fast, normal, x);
    end
end

function res=fwddiff(x)
    res=zeros([size(x), 2]);
    if isa(x, 'cvx')
        res=cvx(res);
    end
    res(1:(end-1), :, 1)=x(2:end, :)-x(1:(end-1), :);
    res(:, 1:(end-1), 2)=x(:, 2:end)-x(:, 1:(end-1));
end

function res=fwddiff_a(x)
    res=zeros(size(x));
    if isa(x, 'cvx')
        res=cvx(res);
    end
    res(1:(end-1), :)=                  x(2:end, :)+x(1:(end-1), :);
    res(:, 1:(end-1))=res(:, 1:(end-1))+x(:, 2:end)+x(:, 1:(end-1));
    res=4*res;
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

function res=bwddiff_a(x)
    res=zeros(size(x));
    if isa(x, 'cvx')
        res=cvx(res);
    end
    res(2:end, :)=              x(2:end, :)+x(1:(end-1), :);
    res(:, 2:end)=res(:, 2:end)+x(:, 2:end)+x(:, 1:(end-1));
    res=4*res;
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

function res=ctrdiff_a(x)
    res=zeros(size(x));
    if isa(x, 'cvx')
        res=cvx(res);
    end

    res(2:(end-1), :)=                  (x(3:end, :)+x(1:(end-2), :))/2;
    res(:, 2:(end-1))=res(:, 2:(end-1))+(x(:, 3:end)+x(:, 1:(end-2)))/2;
    res=4*res;
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
