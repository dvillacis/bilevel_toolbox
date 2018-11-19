function op=make_op_conj(a)
    op.val=@(x) a.conj(x);
    op.conj=@(y) a.val(y);
    op.bound=a.bound;
    op.lowerbound=a.lowerbound;
end
