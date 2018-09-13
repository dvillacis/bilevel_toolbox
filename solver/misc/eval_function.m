function curr_norm = eval_function(fg,Fp,x,s,param)
%EVAL_FUNCTION internal evaluation function
  curr_norm = fg.eval(x);
  for ii = 1:length(Fp)
    curr_norm = curr_norm + Fp{ii}.eval(x);
  end
end
