% Compute the gradient of the huberised one-norm of a vector field.
function hg=hubergrad(g, gamma)
    hg=gamma*g./max(1, gamma*sqrt(sum(g.^2, 3)));
end
