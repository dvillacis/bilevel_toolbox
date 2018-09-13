function algo = algoname(name)

  switch lower(name)
    % case 'forward_backward'
    %   algo = forward_backward_alg();
    case 'chambolle_pock'
      algo = chambolle_pock_alg();
    otherwise
      error('Unknown algorithm name');
  end

end
