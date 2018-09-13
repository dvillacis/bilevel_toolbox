function [sol,info] = chambolle_pock(x_0, f_1, f_2, param)

  param.algo = 'CHAMBOLLE_POCK';
  [sol, info] = solvep(x_0,{f_1, f_2}, param);

end
