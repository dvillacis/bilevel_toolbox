function [grad] = solve_pd_ntr_gradient(u,lambda,dataset,params)
    original = dataset.get_target(1);
    noisy = dataset.get_corrupt(1);

    [M,N] = size(noisy);
    po = PatchOperator(size(lambda),[M,N]);
    
    lambda_out = po.val(lambda);
    nabla = gradient_matrix(M,N);
    A = spdiags(lambda_out(:),0,M*N,M*N); % Build diagonal matrix with parameters to estimate
    B = nabla';
    Ku = nabla*u(:); %Discrete gradient matrix
    nKu = xi(Ku,M,N); %Discrete euclidean norm

    if params.regularized_model == false
        % Get partition active-inactive
        act = (nKu<1e-7); %TODO: Specify a partition of the possible biactive set
        inact = 1-act;
        Act = spdiags(act,0,2*M*N,2*M*N);
        Inact = spdiags(inact,0,2*M*N,2*M*N);

        % Get the adjoint state
        denominador = Inact*nKu+act;
        prodKuKu = outer_product(Ku./(denominador.^3),Ku,M,N);
        C = -Inact*(prodKuKu-spdiags(1./denominador,0,2*M*N,2*M*N))*nabla;
        D = speye(2*M*N);
        E = Act*nabla;
        F = sparse(2*M*N,2*M*N);
        %Adj = [A B;C D;E F];
        Adj = [A B;E-C Inact+sqrt(eps)*Act];
        Track = [u(:)-original(:);sparse(2*M*N,1)];
        mult = Adj\Track;
        adj = mult(1:N*M);
    else
        % Active, weak active and inactive sets \xi(\nabla^h y)>=1.
        gamma = 100;
        act1=gamma*nKu-1;
        act=spones(max(0,act1(1:M*N)-1/(2*gamma)));
        Act=spdiags(act,0,M*N,M*N);
        inact=spones(min(0,act1(1:M*N)+1/(2*gamma)));
        sact=sparse(1-act-inact);
        Sact=spdiags(sact,0,M*N,M*N);

        % Diagonal matrix corresponding to regularization function
        den=(Act+Sact)*nKu(1:M*N)+inact;
        mk=(act+Sact*(1-gamma/2*(1-gamma*nKu(1:M*N)+1/(2*gamma)).^2))./den;
        Dmi=spdiags(kron(ones(2,1),mk+gamma*inact),0,2*M*N,2*M*N);

        % Negative term in the derivative
        subst=spdiags(act+Sact*(1-gamma/2*(1-gamma*nKu(1:M*N)+1/(2*gamma)).^2),0,M*N,M*N);
        subst=kron(speye(2),subst);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Construction of the Hessian components corresponding to each
        % equation in the optimality system
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Matrix of product (Gy)*(Gy)^T for the Newton step divided by
        % the Frobenius norm to the cube

        H4=prodesc(Ku./kron(ones(2,1),den.^2),Ku);

        % Matrix of product (Gy)*(Gy)^T for the Newton step

        prodKuKu=prodesc(Ku,Ku);

        sk2=(Sact*gamma^2*(1-gamma*nKu(1:M*N)+1/(2*gamma)))./(den.^2);
        sk2=spdiags(kron(ones(2,1),sk2),0,2*M*N,2*M*N);

        % Hessian matrix components corresponding to the first equation
        % in the optimality system
        % TODO: This is the second derivative of the l2 norm (check it out!) -> Check derivation!

        hess22=Dmi*nabla-kron(speye(2),(Act+Sact))*Dmi*H4*nabla+sk2*prodKuKu*nabla;

        % Adjoint state is solution of the linear system

        adj=((B*hess22)'+A')\(u(:)-original(:));

    end

    % Calculating the gradient
    beta = 0.1;
    grad = po.conj((noisy-u).*reshape(adj,M,N)) + beta * lambda;
end

function nXi = xi(p,m,n)
  p = reshape(p,m*n,2);
  a = sqrt(sum(p.^2,2));
  nXi = [a;a];
end

function prod = outer_product(p,q,m,n)
  p = reshape(p,m*n,2);
  q = reshape(q,m*n,2);
  a = p(:,1).*q(:,1);
  b = p(:,1).*q(:,2);
  c = p(:,2).*q(:,1);
  d = p(:,2).*q(:,2);
  prod = [spdiags(a,0,m*n,m*n) spdiags(b,0,m*n,m*n); spdiags(c,0,m*n,m*n) spdiags(d,0,m*n,m*n)];
end

%Matriz con el producto q*p^T para el paso de Newton
function P=prodesc(q,p,N)
    if nargin<3
        N=2;
    end
    
    n=size(q,1)/N;

    P=sparse(N*n, N*n);

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


