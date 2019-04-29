% Compute the huberised one-norm of a vector field.
% Parameters:
%    Gy - the vector field, as a 1d vector. 
%    pts - dimension of the vectors of the vector field.
%    gamma - Huber parameter (larger is less regularisation)
function h=huber(Gy, pts, gamma)
    sz=size(Gy, 1);
    N=sz/pts;
    
    % Calculate pointwise Euclidean norm of Gy, expanded
    % to have same dimensions as Gy.
    nGy=xi(Gy, N);

    % We only need first part.
    nGy=nGy(1:pts);
    
    act1=nGy-1/gamma;
    inact=spones(min(0,act1));  %Inactive indicator vector
    act=1-inact;    %Active indicator vector
    
    hh=(gamma/2)*inact.*(nGy.*nGy)+act.*(nGy-1/(2*gamma));
    h=sum(hh);
end
