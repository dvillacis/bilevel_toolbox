function [sol,info] = forward_backward(x_0, f_1, f_2, param)

  param.algo = 'FORWARD_BACKWARD';
  [sol, info] = solvep(x_0,{f_1, f_2}, param);

end
