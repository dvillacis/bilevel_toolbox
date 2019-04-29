% SOLVE_GENERIC_L1_L2 Generic solver for several image processing problems
% This solver receives an abstract initial structure to support different
% image processing models and solves those by using an extension of the
% Hintermüller-Stadler semismooth Newton method; see
%
%   De los Reyes, Schönlieb, and Valkonen, Bilevel parameter learning for
%   higher-order total variation regularisation models,
%   https://arxiv.org/abs/1508.07243, Journal of Mathematical Imaging and
%   Vision 57 (2017), 1–25, doi:10.1007/s10851-016-0662-8
%
% INPUTS
%   lambda: l2 fidelity weights
%   alpha: Huber-l1 fidelity weights
%   Ks: cell array of operators corresponding to l2 terms
%   Bs: cell array of operators corresponding to the l1 terms
%   z: l2 data
%   q: l1 data
%   gamma: Huber regularization parameter for l1 terms
%   param: struct with the algorithm specific parameters
%     .gradstop (boolean, default: false): stop if gradient small
%     .dinit (vector, default: zero): dual variable initailisation
%     .maxiter (integer, default: 100): maximum iteration count
%     .armijo (boolean, default: true): perform Armijo line search
%     .armijo_c (double, default: 1e-4): Armijo line search parameter c.
%     .target_fraction (double, default: 1e-5): step length / norm fractional
%         stoping criterion
%
% OUTPUTS
%   sol: minimizer for the optimization problem
%   val: function values per iteration
%
function [sol,vals] = solve_generic_l1_l2_hs(lambda,alpha,Ks,Bs,z,q,gamma,yinit,param)
    % Default parameters
    param=add_default(param, 'gradstop', false);
    param=add_default(param, 'dinit', []);
    param=add_default(param, 'maxiter', 100);
    param=add_default(param, 'armijo', true);
    param=add_default(param, 'armijo_c', 1e-4);
    param=add_default(param, 'target_fraction', 1e-5);

    % L1 data is not supported at the moment
    for j=1:length(q)
        assert(isempty(q{j}) || norm(q{j}(:))==0);
    end

    % Only square images supported
    assert(size(yinit, 1)==size(yinit, 2));

    n=size(yinit, 1);
    m=n;

    y=yinit(:);
    itercount=0;
    normfrac=1;
    converged=1;
    vals=[];

    % First row of the Hessian.
    GmTs={};
    Gms={};
    for j=1:length(Bs)
        Gms{j}=sparse(Bs{j}.matrix());
        GmTs{j}=transpose(Gms{j});
    end
    B=sparse(cat(2, GmTs{:}));
    HessQ1=0;
    zk=0;
    for k=1:length(Ks)
        thisK=sparse(Ks{k}.matrix());
        HessQ1=HessQ1+(2*lambda{k})*thisK'*thisK;
        % Expand target data into primal domain. Need to vectorise
        zk=zk+(2*lambda{k})*thisK'*(z{k}(:));
    end
    % Number of rows for dual conditions in the Hessian
    drows=size(cat(1, Gms{:}), 1);
    % Calculate first row for dual codition corresponding
    % to each Gms{j} in the Hessian.
    %Gms_rows=cumsum([size(A, 1) cellfun(@(x) size(x, 1), Gms)]);
    % Calculate first column for dual variable corresponding
    % to each Gms{j} in the Hessian.
    Gms_cols=cumsum([1 cellfun(@(x) size(x, 1), Gms)]);

    val=@(yy, dd) calc_val(lambda,alpha,Ks,Bs,z,q,gamma,reshape(yy, n, m),param);
    grad=@(yy, dd) reshape(calc_grad(lambda,alpha,Ks,Bs,z,q,gamma,reshape(yy, n, m),param), size(yy));

    % Initialise dual variable
    if ~isempty(param.dinit)
        d=param.qinit;
    else
        d=zeros(size(B, 2),1);
    end

    Ds={};
    eta2d={};
    Ms={};

    %% Newton loop
    while normfrac>param.target_fraction
        if itercount>param.maxiter
            fprintf(2, 'Reached maximum allowed iterations in denoise_gen (residual=%g).', normfrac);
            converged=0;
            break;
        end
        itercount=itercount+1;

        gg=grad(y, d);
        gn=norm(gg, 2);%/max(1, norm(y, 2));
        if gn<param.target_fraction
            if param.gradstop
                fprintf(1, 'G');
                break;
            else
                fprintf(1, 'g');
            end
        end
        %fprintf(1, '\r%c', iterdisp(1+mod(itercount, length(iterdisp))));
        %fprintf(1, '.');

        % Initialise Hessian and the right-hand-side of
        % the Newton step equation.
        Hess=HessQ1;
        eta2=-(HessQ1*y+B*d)+zk;

        Ds={};
        eta2d={};
        Ms={};
        Hs={};
        e2s={};

        for j=1:length(Gms)
            G=Gms{j};
            a=alpha{j};
            ga=gamma{j};

            %  Evaluate operator G
            Gy=G*y;

            sz=size(Gy, 1);
            pts=n*m;
            N=sz/pts;

            % Calculate pointwise Euclidean norm of Gy, expanded
            % to have same dimensions as Gy.
            nGy=pointwise_norm_replicated(Gy, n, m);

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
            act1=nGy-1/ga;
            inact=spones(min(0,act1));  %Inactive indicator vector
            act=1-inact;    %Active indicator vector
            Act=spdiags(act,0,sz,sz); %Diagonal matrix constructed from active vector

            % Because of max(1, g*nGy) instead of max(1/g, nGy)
            % here, we have the multiplications by g in Hess and eta2
            % below.
            % {Dm here} = {(1/g)Dm for Dm as in the paper}
            % because {our g}={1/g in the paper}.
            mk=max(1, ga*nGy);
            Dm=spdiags(mk,0,sz,sz);

            % Modification of the Hessian: projection of d to a-ball.
            dG=d(Gms_cols(j):Gms_cols(j+1)-1);
            ndG=pointwise_norm_replicated(dG, n, m);
            mqk=max(a,ndG);
            dpG=spdiags(a./mqk,0,sz,sz)*dG;
            Ddp=spdiags(dpG,0,sz,sz);

            % Need to use the fact that Dm is invertible to solve dd,
            % in order to keep the system of reasonable size for TGV.
            %
            % [H      GmTs{1} ...  GmTs{N}] [delta_u ]   [eta2    ]
            % [Ms{1}  Ds{1} .         ] [delta_p1] = [eta2d{1}]
            % [...           .        ] [...     ]   [...     ]
            % [Ms{N}          .  Ds{N}] [delta_pN]   [eta2d{N}]
            %
            % ==> delta_pi=(Ds{i}) (eta2d{i} - Ms{i}delta_u)
            %     H + sum_i GmTs{i}inv(Ds{i}) (eta2d{i} - Ms{i}delta_u) = eta2
            %     ==> (H - sum_i GmTs{i}inv{Ds{i}}Ms{i})delta_u
            %                   = eta2 - sum_i GmTs{i}inv(Ds{i}) eta2d{i}
            %
            % We need the multiplications by g {=1/g in the paper}
            % because of the corresponding scaling of Dm.
            Ds{j}=Dm;
            Ms{j}=ga*Act*Ddp*P*G-ga*a*G;
            eta2d{j}=-Dm*dG+ga*a*Gy;
            Hess=Hess-GmTs{j}*(Dm \ Ms{j});
            eta2=eta2-GmTs{j}*(Dm \ eta2d{j});
        end

        dy=Hess\eta2;
        dd=zeros(size(d));
        for j=1:length(Gms)
            dd(Gms_cols(j):Gms_cols(j+1)-1)=Ds{j} \ (eta2d{j}-Ms{j}*dy);
        end

        v0=val(y, q);
        vals=[vals, v0];

        if ~param.armijo
            normfrac=norm(dy,2)/max(1, norm(y,2));
            y=y+dy;
            d=d+dd;
            st='.';
        else
            normfracbase=norm(dy,2)/max(1, norm(y,2));
            sigma=1;
            g0=dy'*gg;
            maxsigmai=32;
            sigmai=1;
            ytry=y;
            dtry=d;
            while sigma*normfracbase > param.target_fraction
                ytry=y+sigma*dy;
                dtry=d+sigma*dd;
                v=val(ytry, dtry);
                vals(end)=v;
                if v <= v0 + sigma*param.armijo_c*g0
                    break
                end
                sigma=sigma/2;
                sigmai=sigmai+1;
            end
            y=ytry;
            d=dtry;
            normfrac=sigma*normfracbase;
            if sigmai==1
                st='.';
            elseif sigmai>=maxsigmai
                st='!';
            else
                st='_';
            end
        end

        fprintf(1, st);
    end

    %fprintf(1, '\r \r');
    sol=reshape(y, [n, m]);
end

function val=calc_val(lambda,alpha,Ks,Bs,z,q,gamma,y,param);
    val=0;
    for j=1:length(Ks)
        val=val+(2*lambda{j})*norm2sq(Ks{j}.val(y)-z{j})/2;
    end
    for j=1:length(Bs)
        val=val+(2*alpha{j})*huber(Bs{j}.val(y)-q{j}, gamma{j});
    end
end

function grad=calc_grad(lambda,alpha,Ks,Bs,z,q,gamma,y,param)
    grad=0;
    for j=1:length(Ks)
        grad=grad+lambda{j}*Ks{j}.conj((Ks{j}.val(y)-z{j}));
    end
    for j=1:length(Bs)
        grad=grad+alpha{j}*Bs{j}.conj(hubergrad(Bs{j}.val(y), gamma{j}));
    end
end


