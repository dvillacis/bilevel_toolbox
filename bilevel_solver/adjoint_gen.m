%
% Calculate the solution p to the adjoint equation, as well 
% as the gradient of the target functional g in terms of u 
% through the adjoint equation.
%
% The parameters are:
%   A - Weighted Laplacian for the Hilbert space diffusion term.
%   K - Matrix for extracting the interesting part from the expanded image (used to handle TGV etc.
%		in a uniform manner).
%   Bs - Cell array of divergence (negative?) for the sum of regularisation terms; see the paper.
%   Gs - Corresponding gradients (TODO: only one of these should be needed!)
%   n - dimension of the assumed-square image.
%   u_fid - weight of fidelity term.
%   u_reg - weight vector for the regularisation terms.
%   gammas - vector of huber gammas.
%   z - ground-truth image
%   yn - noisy image
%   y - current iterate (expanded image containing extra variables
%       for TGV2 etc.; K*y is the image itself).
%   FIDSOL - 1 of fidelity weight is unknown, zero otherwise
%   REGSOL - vector of ones and zeros corresponding to unknown
%            regularisation term weights
%   eta - adjoint state(?) (negative cost functional gradient at current weight)

function [p, gp]=adjoint_gen(A,K,Bs,Gs,n,u_fid,u_reg,gammas,z,yn,y,FIDSOL,REGSOL, eta)
    m=n;
    g=1;

    hess=A;
    %greg=zeros(size(y, 1), length(Gs));
    hgammag={};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fidelity part
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % lambda_i phi_i''(u)
    hess=hess+u_fid*K'*K;
    %eta=-K'*(K*y-z);
    %gfid=K'*(K*y-yn);
    gfidx=(K*y-yn);
    
    hs={};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Regulariser part
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %parfor i=1:length(Gs)
    for i=1:length(Gs)
        G=Gs{i};
        B=Bs{i};
        gama=gammas(i);
        c=gama;
        
        Gy=G*y;
        sz=size(Gy, 1);
        pts=m^2;
        N=sz/pts;

        % Calculate pointwise Euclidean norm of Gy, expanded
        % to have same dimensions as Gy.
        nGy=xi(Gy, N);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Divergence rowwise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % This function computes the discrete euclidian norm of gradient of y, \xi(\nabla^h y). 

        % Active, weak active and inactive sets \xi(\nabla^h y)>=g.
        act1=gama*nGy-g;
        act=spones(max(0,act1(1:m^2)-1/(2*c)));
        Act=spdiags(act,0,m^2,m^2);
        inact=spones(min(0,act1(1:m^2)+1/(2*c)));
        Inact=spdiags(inact,0,m^2,m^2);
        sact=sparse(1-act-inact);
        Sact=spdiags(sact,0,m^2,m^2);

        %Matriz diagonal correspondiente a la
        %funcion de regularizacion
        denominador=(Act+Sact)*nGy(1:m^2)+inact;
        mk=(g*act+Sact*(g-c/2*(g-gama*nGy(1:m^2)+1/(2*c)).^2))./denominador;
        Dmi=spdiags(kron(ones(N,1),mk+gama*inact),0,sz,sz);

        %Termino negativo en la derivada
        subst=spdiags(g*act+Sact*(g-c/2*(g-gama*nGy(1:m^2)+1/(2*c)).^2),0,m^2,m^2);
        subst=kron(speye(N),subst);

        %Matriz para condicionar la Hessiana
        mk=max(g,gama*nGy);
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

        sk2=(Sact*gama^2*(g-gama*nGy(1:m^2)+1/(2*c)))./(denominador.^2);
        sk2=spdiags(kron(ones(N,1),sk2),0,sz,sz); 
                   
        % max(g,gamma*N(E^h y)). 
        mk=max(g,gama*nGy);
        Dm=spdiags(mk,0,sz,sz);

                
        %Hess1: Componentes de la matriz hessiana correspondientes a la
        % primera ecuacion en el sistema de optimalidad
        hess22=Dmi-kron(speye(N),(Act+Sact))*Dmi*h41baux+sk2*prodGyGy;
        
        % G^* h_gamma'(Gy) G
        hs{i}=u_reg(i)*B*hess22*G;
        % G^* h_gamma(G y)
        %greg(:, i)=B*Dmi*Gy;
        hgammag{i}=Dmi*Gy;
    end

    for i=1:length(Gs)
        hess=hess+hs{i};
    end
    
    p=hess \ eta;
    
    gp=zeros(sum([FIDSOL; REGSOL]), 1);
    idx=1;
    if FIDSOL
        %gp(idx)=multi_iprod(gfid, p, m);
        gp(idx)=multi_iprod(gfidx, K*p, m);
        idx=idx+1;
    end
    for i=1:length(Gs)
        if REGSOL(i)
            %gp(idx)=multi_iprod(greg(:,i), p, m);
            gp(idx)=multi_iprod(hgammag{i}, Gs{i}*p, m);
            idx=idx+1;
        end
    end
end
