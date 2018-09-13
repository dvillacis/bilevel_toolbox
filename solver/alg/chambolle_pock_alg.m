function s = chambolle_pock_alg()
  s.name = 'CHAMBOLLE_POCK';
  s.initialize = @(x_0, fg, Fp, param) chambolle_pock_initialize(x_0,fg,Fp,param);
  s.algorithm = @(x_0,fg,Fp,sol,s,param) chambolle_pock_algorithm(Fp,sol,s,param);
  s.finalize = @(x_0,fg,Fp,sol,s,param) sol;
end

function [sol,s,param] = chambolle_pock_initialize(x_0,fg,Fp,param)

  if numel(Fp) == 1
    Fp{2}.prox = @(x,T) x;
    Fp{2}.eval = @(x) eps;
  end

  s = struct;
  isL = 1;
  if isfield(Fp{1},'L')
    s.ind = [2,1];
    L = Fp{1}.L;
    Lt = Fp{1}.Lt;
  elseif isfield(Fp{2},'L')
    s.ind = [1,2];
    L = Fp{2}.L;
    Lt = Fp{2}.Lt;
  else
    L = @(x) x;
    Lt = @(x) x;
    s.ind = [1,2];
    isL = 0;
  end

  if isfield(Fp{s.ind(2)}, 'norm_L')
    s.norm_L = Fp{s.ind(2)}.norm_L;
  else
    s.norm_L = 1;
    if isL
      warning('You should give f.norm_L = ||L||^2. Setting it to 1!');
    end
  end

  % test for timesteps
  if ~isfield(param,'sigma') && ~isfield(param,'tau')
    s.tau = 1/sqrt(s.norm_L);
    s.sigma = 1/sqrt(s.norm_L);
  elseif isfield(param,'tau')
    s.tau = param.tau;
    s.sigma = 1/(s.tau * s.norm_L);
    if sigma < 0
      error('Tau is too big!');
    end
  elseif isfield(param,'sigma')
    s.sigma = param.sigma;
    s.tau = 1/(s.sigma * s.norm_L);
  else
    s.sigma = param.sigma;
    s.tau = param.tau;
  end

  s.OpL = L;
  s.OpLt = Lt;
  s.dual_var = s.OpL(x_0);

  sol = x_0;
  s.x_n = x_0;

  if fg.beta
    error('CHAMBOLLE POCK needs only function with proximal operators');
  end

  if param.lambda >= 1
    warning('param.lambda is greater than 1!')
  end
end

function [sol,s] = chambolle_pock_algorithm(Fp,sol,s,param)

  if numel(Fp)==1
    Fp{2}.prox = @(x,T) x;
    Fp{2}.eval = @(x) eps;
  end

  % Algorithm
  s.dual_var = prox_adjoint( s.dual_var + s.sigma * s.OpL(sol), s.sigma, Fp{s.ind(2)});
  x_n_old = s.x_n;
  s.x_n = Fp{s.ind(1)}.prox(s.x_n - s.tau * s.OpLt(s.dual_var), s.tau);
  sol = s.x_n + param.lambda * (s.x_n - x_n_old);
end
