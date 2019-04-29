% Compute the huberised one-norm of a vector field.
function h=huber(g, gamma)
    ngsq=sum(g.^2, 3);
    ng=sqrt(ngsq);

    act1=ng-1/gamma;
    inact=spones(min(0,act1));  %Inactive indicator vector
    act=1-inact;    %Active indicator vector

    hh=(gamma/2)*inact.*(ngsq)+act.*(ng-1/(2*gamma));
    h=sum(hh(:));
end
