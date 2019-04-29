% Compute the gradient of the huberised one-norm of a vector field.
% Parameters:
%    Gy - the vector field, as a 1d vector. 
%    pts - dimension of the vectors of the vector field.
%    gamma - Huber parameter (larger is less regularisation)
function hg=hubergrad(Gy, pts, gamma)
    sz=size(Gy, 1);
    N=sz/pts;
    
    nGy=xi(Gy, N);
    
    hg=gamma*Gy./max(1, gamma*nGy);
end
