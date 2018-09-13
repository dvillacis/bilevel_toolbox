function solver = select_solver(Fgrad, Fprox)

  n_linear = numberofL(Fprox);
  n_prox = numel(Fprox);
  n_grad = numel(Fgrad);

  if n_grad > 0 % There are smooth functions
    if n_linear > 2
      error('Sorry, no solver is able to solve your problem yet!');
    end
    if n_prox == 0
      solver = 'GRADIENT_DESCENT';
    elseif n_prox == 2
      solver = 'FB_BASED_PRIMAL_DUAL';
    elseif (n_prox <= 2) && (n_linear > 0)
      solver = 'FB_BASED_PRIMAL_DUAL';
    elseif n_prox == 1
      solver = 'FORWARD_BACKWARD';
    else
      solver = 'GENERALIZED_FORWARD_BACKWARD';
    end
  else % All functions are nonsmooth
    if (n_prox <= 2) && (n_linear == 1)
      solver = 'CHAMBOLLE_POCK';
    elseif (n_linear > 1)
      solver = 'SDMM';
    elseif n_prox == 2
      solver = 'DOUGLAS_RACHFORD';
    else
      solver = 'PPXA';
    end
  end
end

function n = numberofL(Fp)
  %NUMBEROFL Returns the number of functions with a linear operator inside
  n = 0;
  for ii = 1:length(Fp)
    if isfield(Fp{ii},'L')
      n = n + 1;
      if ~isfield(Fp{ii},'Lt')
        warning('You did not define the Lt operator!');
      end
    end
  end
end
