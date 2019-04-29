%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEMISMOOTH NEWTON METHOD FOR THE SOLUTION OF A GENERAL
% DENOISING MODEL, DISCRETIZATION BY FINITE DIFFERENCES,
% SIMILAR TO TV CASE IN (Hintermueller and Stadler (2004))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y1, converged, cont, q1]=denoise_gen(A,K,Bs,Gs,n,g,alphas,gammas,noise,yinit, qinit, gradstop)
    %% Notation:
    %   A: corresponds to the stiffness matrix
    %   K: corresponds to K matrix
    %   Bs: cell array of matrices, conjugate to Gs
    %   G: cell array of matrices, corresponding to multiple gradient operators
    %   n: is the numer of disc. points
    %   g: fidelity weight
    %   alphas: vector of regularisation weights
    %   mu: is the epsilon parameter in front of the Laplacian
    %   gamma: is the Huber regularization parameter
    %   noise: is the noisy image
    %   yinit: the the image from which the method starts

    %% Initialization
    %g2=g.*ones(n^2,1);          %We built a vector from the parameter
    %U=spdiags(g2,0,n^2,n^2);    %Diagonal matrix with coefficient in the diagonal
    %g=1;
    %c=gamma;
    m=n;
    y=yinit;
    cont=0;
    norma=1;
    rate=0;
    residuo=0;
    maxiter=100;
    converged=1;
    armijo=1;
    armijo_c=1e-4;
    tgt_norma=1e-5;
    %norq=0;
    if nargin<12
        gradstop=false;
    end


    % First row of the Hessian.
    B=sparse(cat(2, Bs{:}));
    HessQ1=[A+g*sparse(K'*K)];
    % Number of rows for dual conditions in the Hessian
    drows=size(cat(1, Gs{:}), 1);
    % Calculate first row for dual codition corresponding 
    % to each Gs{j} in the Hessian.
    %Gs_rows=cumsum([size(A, 1) cellfun(@(x) size(x, 1), Gs)]);
    % Calculate first column for dual variable corresponding 
    % to each Gs{j} in the Hessian.
    Gs_cols=cumsum([size(A, 2) cellfun(@(x) size(x, 1), Gs)])-length(y)+1;
    cols_p=[1:size(A, 2)];

    val=@(yy, qq) calc_value(K, A, Gs, n, g, alphas, gammas, Gs_cols, noise, yy, qq);
    grad=@(yy, qq) calc_grad(K, A, Gs, n, g, alphas, gammas, Gs_cols, noise, yy, qq);
    
    % Expand target data into primal domain
    z=g*K'*noise;

    % Initialise dual varibale
    if ~isempty(qinit)
        q=qinit;
    else
        q=zeros(size(B, 2),1);
    end
    
    Ds={};
    eta2q={};
    Ms={};
    
    %iterdisp='/-\|';
    
    %% Newton loop
    while norma>tgt_norma
        if cont>maxiter
            fprintf(2, 'Reached maximum allowed iterations in denoise_gen (residual=%g).', norma);
            converged=0;
            break;
        end
        cont=cont+1;
        
        gg=grad(y, q);
        gn=norm(gg, 2);%/max(1, norm(y, 2));
        if gn<tgt_norma
            if gradstop
                fprintf(1, 'G');
                break;
            else
                fprintf(1, 'g');
            end
        end
        %fprintf(1, '\r%c', iterdisp(1+mod(cont, length(iterdisp))));
        %fprintf(1, '.');
        
        % Initialise Hessian and the right-hand-side of
        % the Newton step equation.
        %Hess=[HessR1; sparse(drows, size(HessR1, 2))];
        %eta2=[-HessR1*[y; q]+z; zeros(drows, 1)];
        Hess=HessQ1;
        eta2=-(HessQ1*y+B*q)+z;
        %eta2=-HessQ1*y+z;
        
        Ds={};
        eta2q={};
        Ms={};
        Hs={};
        e2s={};

        %parfor j=1:length(Gs)
        for j=1:length(Gs)
            G=Gs{j};
            alpha=alphas(j);
            gamma=gammas(j);
            
            %  Evaluate operator G
            Gy=(G*y);
            
            sz=size(Gy, 1);
            pts=m^2;
            N=sz/pts;
            
            % Calculate pointwise Euclidean norm of Gy, expanded
            % to have same dimensions as Gy.
            nGy=xi(Gy, N); 
            
            % Here we calculate the generalized differential P of \xi.
            GyD=sparse(pts, N*pts);
            for i=1:N
                idx=1+(i-1)*pts:i*pts;
                GyD(:, idx)=spdiags(Gy(idx),0,pts,pts);
            end
            R=kron(ones(N,1), GyD);
            Nin=spdiags(1./nGy,0,sz,sz);
            P=Nin*R;
            
            % Active and inactive components \xi(\nabla^h y)>=g.
            %
            act1=nGy-1/gamma;
            inact=spones(min(0,act1));  %Inactive indicator vector
            act=1-inact;    %Active indicator vector
            Act=spdiags(act,0,sz,sz); %Diagonal matrix constructed from active vector
            
            % Because of max(1, gamma*nGy) instead of max(1/gamma, nGy)
            % here, we have the multiplications by gamma in Hess and eta2
            % below.
            % {Dm here} = {(1/gamma)Dm for Dm as in the paper}
            % because {our gamma}={1/gamma in the paper}.
            mk=max(1, gamma*nGy); 
            Dm=spdiags(mk,0,sz,sz); 
            
            % Modification of the Hessian: projection of q to alpha-ball.
            qG=q(Gs_cols(j):Gs_cols(j+1)-1);
            nqG=xi(qG, N);
            mqk=max(alpha,nqG);
            qpG=spdiags(alpha./mqk,0,sz,sz)*qG;
            Dqp=spdiags(qpG,0,sz,sz);
            
            %% Construct further line of linear system for 
            %% regulariser |Gs{j}y|
            %Hess(rows_G, cols_p)=gamma*Act*Dqp*P*G-gamma*alpha*G;
            %Hess(rows_G, cols_G_d)=Dm;
            %eta2(rows_G)=[-Dm*qG+gamma*alpha*Gy];
            
            % Need to use the fact that Dm is invertible to solve dq, 
            % in order to keep the system of reasonable size for TGV.
            %
            % [H      Bs{1} ...  Bs{N}] [delta_u ]   [eta2    ]
            % [Ms{1}  Ds{1} .         ] [delta_p1] = [eta2q{1}]
            % [...           .        ] [...     ]   [...     ]
            % [Ms{N}          .  Ds{N}] [delta_pN]   [eta2q{N}]
            %
            % ==> delta_pi=(Ds{i}) (eta2q{i} - Ms{i}delta_u)
            %     H + sum_i Bs{i}inv(Ds{i}) (eta2q{i} - Ms{i}delta_u) = eta2
            %     ==> (H - sum_i Bs{i}inv{Ds{i}}Ms{i})delta_u 
            %                   = eta2 - sum_i Bs{i}inv(Ds{i}) eta2q{i}
            %
            % We need the multiplications by gamma {=1/gamma in the paper}
            % because of the corresponding scaling of Dm.
            Ds{j}=Dm;
            Ms{j}=gamma*Act*Dqp*P*G-gamma*alpha*G;
            eta2q{j}=-Dm*qG+gamma*alpha*Gy;
            %Hs{j}=-Bs{j}*(Dm \ Ms{j});
            %e2s{j}=-Bs{j}*(Dm \ eta2q{j});
            Hess=Hess-Bs{j}*(Dm \ Ms{j});
            eta2=eta2-Bs{j}*(Dm \ eta2q{j});
        end
        
        % Solve primal and dual step
        %dx=Hess\eta2;
        %dy=dx(1:length(y));
        %dq=dx(length(y)+1:end);
        
        %for j=1:length(Gs)
        %    Hess=Hess+Hs{j};
        %    eta2=eta2+e2s{j};
        %end
        
        dy=Hess\eta2;
        %dq=blkdiag(Ds{:}) \ (cat(1, eta2q{:})-cat(1, Ms{:})*dy);
        dq=zeros(size(q));
        for j=1:length(Gs)
            dq(Gs_cols(j):Gs_cols(j+1)-1)=Ds{j} \ (eta2q{j}-Ms{j}*dy);
        end

        % Add Armijo line search here?

        
        if ~armijo
            norma=norm(dy,2)/max(1, norm(y,2)); 
            y=y+dy;
            q=q+dq;
            st='.';
        else
            normabase=norm(dy,2)/max(1, norm(y,2)); 
            sigma=1;
            v0=val(y, q);
            g0=dy'*gg;
            maxsigmai=32;
            sigmai=1;
            ytry=y;
            qtry=q;
            while sigma*normabase > tgt_norma
                ytry=y+sigma*dy;
                qtry=q+sigma*dq;
                v=val(ytry, qtry);
                if v <= v0 + sigma*armijo_c*g0
                    break
                end
                sigma=sigma/2;
                sigmai=sigmai+1;
            end
            %fprintf(1, '[sigma=%g;%g,%g,%g]', sigma, v0, v, val(y+dy, q+dq));%armijo_c*g0);
            y=ytry;
            q=qtry;
            norma=sigma*normabase;
            if sigmai==1
                st='.';
            elseif sigmai>=maxsigmai
                st='!';
            else
                st='_';
            end
        end
        
        fprintf(1, st);
        %norma=norm(dy,2)/max(1, norm(y,2)); 
    end

    %fprintf(1, '\r \r');
    
    y1=y;
    q1=q;
end


function val=calc_value(K, A, Gs, m, g, alphas, gammas, Gs_cols, noise, y, q)
    l2normUSC=@(x) sqrt(sum(x.^2));
    val=(g/2)*l2normUSC(K*y-noise)^2+(1/2)*y'*A*y;
    for j=1:length(Gs)
        G=Gs{j};
        val=val+alphas(j)*huber(G*y, m^2, gammas(j));
    end
    if 0
    % Dual; cheating.
    val=val+l2normUSC(noise)^2/2;
    Ksn=K'*noise;
    for j=1:length(Gs)
        G=Gs{j};
        qG=q(Gs_cols(j):Gs_cols(j+1)-1);
        %KsnG=Ksn(Gs_rows(j):Gs_rows(j+1)-1);
        val=val-gammas(j)/(2*alphas(j))*l2normUSC(qG)^2;
        % Should be H^{-1} norm!!
        size(Gs{j}'*qG)
        size(Ksn)
        val=val-(1/2)*l2normUSC(Gs{j}'*qG-Ksn)^2;
    end
    end
end

function grad=calc_grad(K, A, Gs, m, g, alphas, gammas, Gs_cols, noise, y, q)
    grad=g*K'*(K*y-noise)+A*y;
    for j=1:length(Gs)
        G=Gs{j};
        grad=grad+alphas(j)*G'*hubergrad(G*y, m^2, gammas(j));
    end
end

%% Graphic interface (if necessary)
%  figure(1);
%  subplot(2,2,1)
%  y2=reshape(full(y),m,m);
%  surf(y2);
%  view(2);
%  title('Denoised image');
%  axis off;
%  %hold off;
% % 
%  subplot(2,2,2)
%  w2=reshape(full(noise),m,m);
%  surf(w2);
%  view(2);
%  title('Noised image');
%  axis off;
%  
%  subplot(2,2,3)
%  load imagingdata.mat o*
%  orig=original(1:m,1:m);
%  surf(reshape(full(orig),m,m));
%  view(2);
%  title('Original image');
%  axis off;
% % hold off;

 
% CONVERGENCE: residuum, convergence rate and size of active set
% figure(2);
% subplot(2,2,1)
% plot(rate(5:cont+1));
% title('Convergence rate');
% 
% subplot(2,2,2)
% semilogy(residuo);
% title('Residuum');
% 
% subplot(2,2,3)
% plot(siac);
% title('Size of Active set vs. iteration');
