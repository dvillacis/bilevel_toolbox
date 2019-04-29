%
% Calculate the solution p to the adjoint equation, as well as the gradient
% of the target functional g in terms of the parameters alpha and lambda
% (only for those components marked as unknowns through FIDSOL and REGSOL).
%
% The parameters are:
%   Ks: cell array of operators corresponding to l2 terms
%   Bs: cell array of operators corresponding to the l1 terms
%   Gs - Corresponding gradients (TODO: only one of these should be needed!)
%   n - dimension of the assumed-square image.
%   lambda - weight of fidelity term.
%   alpha - weight vector for the regularisation terms.
%   gamma - vector of Huber gammas.
%   f - ground-truth images
%   z - noisy images
%   y - current iterate
%   FIDSOL - vector of ones and zeros corresponding to unknown
%            fidelity term weights
%   REGSOL - vector of ones and zeros corresponding to unknown
%            regularisation term weights
%   eta - adjoint state
function [p, glambda, galpha]=adjoint_gen(Ks,Bs,lambda,alpha,gamma,f,z,FIDSOL,REGSOL,eta)
    n=size(y, 1);
    m=n;

    % Only square images supported at the moment
    assert(m==n);

    g=1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fidelity part
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hess=0;
    gfid={};
    for j=1:length(Ks)
        Km=Ks{j}.matrix();
        hess=hess+(lambda{j}*2)*Km'*Km;
        gfid{j}=(Ks{j}.val(y)-z{j});
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Regulariser part
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    greg={};
    hs={};
    for i=1:length(Gs)
        Gm=Bs{i}.matrix();
        ga=gamma{i};
        c=ga;

        Gy=Gm*y(:);
        sz=size(Gy, 1);
        pts=m*n;
        N=sz/pts;

        % Calculate pointwise Euclidean norm of Gy, expanded to have same dimensions as Gy.
        nGy=pointwise_norm_replicated(Gy, N);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Divergence rowwise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Active, weak active and inactive sets \xi(\nabla^h y)>=g.
        act1=ga*nGy-g;
        act=spones(max(0,act1(1:sz)-1/(2*c)));
        Act=spdiags(act,0,sz,sz);
        inact=spones(min(0,act1(1:sz)+1/(2*c)));
        Inact=spdiags(inact,0,sz,sz);
        sact=sparse(1-act-inact);
        Sact=spdiags(sact,0,sz,sz);

        %Matriz diagonal correspondiente a la
        %funcion de regularizacion
        denominador=(Act+Sact)*nGy(1:sz)+inact;
        mk=(g*act+Sact*(g-c/2*(g-ga*nGy(1:sz)+1/(2*c)).^2))./denominador;
        Dmi=spdiags(kron(ones(N,1),mk+ga*inact),0,sz,sz);

        %Termino negativo en la derivada
        subst=spdiags(g*act+Sact*(g-c/2*(g-ga*nGy(1:sz)+1/(2*c)).^2),0,sz,sz);
        subst=kron(speye(N),subst);

        %Matriz para condicionar la Hessiana
        mk=max(g,ga*nGy);
        Dm=spdiags(mk,0,sz,sz);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%Construction of the Hessian components corresponding to each
        %%equation in the optimality system
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Matriz correspondiente al producto (Gy)*(Gy)^T para el paso de
        %Newton dividida por la norma de frobenius al cubo
        h41baux=prodesc(Gy./kron(ones(N,1), denominador.^2), Gy, N);

        %Matriz correspondiente al producto (Gy)*(Gy)^T para el paso de
        %Newton
        prodGyGy=prodesc(Gy, Gy, N);

        sk2=(Sact*ga^2*(g-ga*nGy(1:sz)+1/(2*c)))./(denominador.^2);
        sk2=spdiags(kron(ones(N,1),sk2),0,sz,sz);
        mk=max(g,ga*nGy);
        Dm=spdiags(mk,0,sz,sz);

        %Hess1: Componentes de la matriz hessiana correspondientes a la
        % primera ecuacion en el sistema de optimalidad
        hess22=Dmi-kron(speye(N),(Act+Sact))*Dmi*h41baux+sk2*prodGyGy;

        hess=hess+alpha{i}*Gm'*hess22*Gm;
        greg{i}=reshape(Dmi*Gy, [n, m, N];
    end

    p = reshape(hess \ reshape(eta, [n*m, 1]), size(y));

    glambda={};
    for j=1:length(Ks)
        if FIDSOL(j)
            glambda{j}=iprod(gfid{j}, Ks{j}.val(p));
        end
    end
    galpha={};
    for i=1:length(Gs)
        if REGSOL(i)
            galpha{i}=iprod(greg{i}, Gs{i}.val(p));
        end
    end
end

function t=iprod(f, p)
    t=sum(f(:).*p(:));
end

%Matriz con el producto q*p^T para el paso de Newton
function P=prodesc(q,p,N)
    n=size(q,1)/N;

    P=sparse(size(q, 1), size(q, 1));

    for j=1:N
        qidx=1+(j-1)*n:j*n;
        qj=q(qidx, 1);
        for i=1:N
            pidx=1+(i-1)*n:i*n;
            pi=p(pidx, 1);

            P(qidx, pidx)=spdiags(pi.*qj, 0, n, n);
        end
    end
end

