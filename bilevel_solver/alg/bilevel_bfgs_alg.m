function s = bilevel_bfgs_alg()
  s.name = 'BFGS';
  s.initialize = @(initial_parameter, lower_level_problem, upper_level_problem, param) bilevel_bfgs_initialize(initial_parameter,lower_level_problem,upper_level_problem,param);
  s.algorithm = @(initial_parameter,lower_level_problem,upper_level_problem,sol,s,param) bilevel_bfgs_algorithm(lower_level_problem,upper_level_problem,sol,s,param);
  s.finalize = @(initial_parameter,lower_level_problem,upper_level_problem,sol,s,param) sol;
end

function [sol,state,param] = bilevel_bfgs_initialize(initial_parameter,lower_level_problem,upper_level_problem,param)
  sol = initial_parameter;
  state.sol_history = sol;
  state.res = 1;
  state.u_history = lower_level_problem.solve(sol,upper_level_problem.dataset);
  state.BFGS = eye(size(state.u_history(1),1)); %Initial BFGS Matrix

  % Default param values
  if ~isfield(param, 'armijo_c')
      param.armijo_c = 1e-4;
  end

  if ~isfield(param, 'wolfe_c')
      param.wolfe_c = 0.1;
  end

end

function [sol,state] = bilevel_bfgs_algorithm(lower_level_problem,upper_level_problem,sol,state,param)

    % Solving the state equation (lower level solver)
    u = lower_level_problem.solve(sol,upper_level_problem.dataset);
    state.u_history = cat(3,state.u_history,u);

    % Getting current cost
    cost = upper_level_problem.eval(u,sol,upper_level_problem.dataset);

    % Saving cost history
    if ~isfield(s, 'l2_cost_history')
        state.l2_cost_history = cost;
    else
        state.l2_cost_history = [s.l2_cost_history cost];
    end

    % Solving the gradient
    state.grad = upper_level_problem.gradient(u,sol,upper_level_problem.dataset,gradient_parameters);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Construct BFGS matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if use_bfgs && iter >= 2
        s = sol -state.solprev;
        t=state.BFGS*s;
        r=(gp+weightregfngrad(u))-(gpprev+weightregfngrad(uprev));
        if (s'*r) < 0
            fprintf(2, '<<NEGATIVE CURVATURE - Skipping BFGS update>>')
        else
            state.BFGS=state.BFGS-1/(s'*t)*kron(t',t)+1/(s'*r)*kron(r',r);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve search direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if use_bfgs
        if find(isnan(BFGS))
            warning('BFGS matrix became NaN');
            error='BFGS matrix';
            erriter=iter;
            return;
        end
        delta=-BFGS\(gp+weightregfngrad(u));
    else
        delta=-(gp+weightregfngrad(u));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Store previous iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    state.solprev = sol;
    state.gradprev = state.grad;
    state.uprev = u;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Line search & state equation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigma=1;
    gg=delta'*(gp+weightregfngrad(uprev));

    mxx=-(uprev-ua)./delta;
    sigma=min([sigma; mxx(mxx>0)/2]);
    mxx=(ub-uprev)./delta;
    sigma=min([sigma; mxx(mxx>0)/2]);

    u_ls0=u;
    y_ls0=y;
    q_ls0=q;
    value_ls0=value;
    res_ls0=res;

    lsiter=0;
    while true
        %u=max(ua,min(uprev+sigma*delta,ub));
        u=uprev+sigma*delta;
        if iscell(y)
            value=weightregfn(u);
            for i=1:length(y)
                fprintf(1, '(IM%d/%d)', i, length(y));
                [y_, dnconverged_,its_, q_]=denoise_gen(A,K,Bs,Gs,m,FIDW(u),REGW(u),gamma,yn{i},y{i},q{i});
                y{i}=y_;
                q{i}=q_;
                dnconverged(i)=dnconverged_;
                its(i)=its_;
                value=value+costfn(y_, z{i});
            end
        else
            [y, dnconverged,its, q]=denoise_gen(A,K,Bs,Gs,m,FIDW(u),REGW(u),gamma,yn,y,q);
            value=costfn(y, z)+weightregfn(u);
        end

        if ~isempty(find(~dnconverged)) && isempty(error)
            error='denoise did not converge';
            erriter=iter;
        end

        dniters=dniters+its;
        res=norm((u-uprev)./u);%+max(0, valprev-value)/valprev;

        if value<value_ls0
            u_ls0=u;
            y_ls0=y;
            q_ls0=q;
            value_ls0=value;
            res_ls0=res;
        end

        if res<tgt_res
            % Too small sigma.
            if value_ls0 < valprev %iter < 10
                % The BFGS matrix B may still be bad - give a chance
                %warning(' line search reached too small residual!!! - reverting to best step');
                fprintf(1, 'R');
                u=u_ls0;
                y=y_ls0;
                q=q_ls0;
                value=value_ls0;
                res=res_ls0;
            else
                %warning(' line search reached too small residual!!! - giving up');
                fprintf(1, 'L');
                y=yprev;
                q=qprev;
                p=pprev;
                value=valprev;
                res=0;
                %if isempty(error)
                %    error='line search residual small';
                %    erriter=iter;
                %end
            end
            break
        end

        switch Armijo
            case 0 % No line search
                break
            case 1 % Armijo
                if value <= valprev + sigma*armijo_c*gg
                    break
                end
                %uprev-u
                %gg
                %[value, valprev, valprev + sigma*armijo_c*gg]
                sigma=sigma/2;
            case 2 % Wolfe
                if iscell(y)
                    error('Wolfe unimplemented for multiple images')
                end
                [pN, gpN]=adjoint_gen(A,K,Bs,Gs,m,FIDW(u),REGW(u),gamma,z,yn,y,FIDSOL,REGSOL, eta);
                ggN=delta'*(gpN+weightregfngrad(u));
                %fprintf(1, '[gg=%f, ggN=%f]', gg, ggN);
                if value <= valprev + sigma*armijo_c*gg && ggN >= wolfe_c*gg
                    %fprintf(1, 'W');
                    break
                end
                sigma=sigma/2;
            case 3 % Quadratic
                if value <= valprev + sigma*armijo_c*gg
                    break
                end
                sigma=-gg*sigma^2/(2*(value-valprev-gg*sigma));
                %sigma=min(sigmaq, sigma/2);
            case 4 % phi(sigma)=c1-c2*sigma*exp(-c3*sigma)
                if value <= valprev + sigma*armijo_c*gg
                    break
                end
                c1=valprev;
                c2=-abs(gg);
                c3=(1/sigma)*(log(sigma)-log(abs((c1-value)/c2)));
                if c3<0
                    sigma=sigma/2;
                else
                    sigma=min(1/c3, sigma/2);
                end
        end

        ustr=num2str(u', '%f, ');
        fprintf(1, '<val=%f,u=%s>', value, ustr);
        fprintf(1, '/')
        lsiter=lsiter+1;
    end
    fprintf(1, '\n');

    vals=[vals, value];
    residuo=[residuo,res];

end
