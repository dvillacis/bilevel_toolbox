% NONSMOOTH_TRUST_REGION_ALG Trust-Region Algorithm for nonconvex-nonsmooth
% functions
% This solver receives an abstract lower and upper level problem descriptions
% of a Bilevel problem and finds a Clarke stationary point by means of
% a two stages trust-region model; see
%
%   De los Reyes, Villacis, Bilevel Parameter Learning for Image Denoising,
%   https://arxiv.org/abs/XXXX,
%
% INPUTS
%   x_0: Initial parameter
%   lower_level_problem: Lower Level Problem struct instance
%   upper_level_problem: Upper Leve� Problem struct instance
%   param: struct with the algorithm specific parameters
%     .minradius (real, default: 1e-4): model switching radius
%     .use_bfgs (boolean, default: false): use a second order model using
%                                           bfgs approximation
%     .use_lbfgs (boolean, default: false): use a second order model using
%                                           a limited memory bfgs
%     .use_linear (boolean, default: true): use a linear model
%
% OUTPUTS
%   sol: current state of the parameter solution
%   info: algorithm run information


function s = nonsmooth_trust_region_alg()
    s.name = 'NONSMOOTH_TRUST_REGION';
    s.initialize = @(x_0, lower_level_problem, upper_level_problem, param) nonsmooth_trust_region_initialize(x_0,lower_level_problem,upper_level_problem,param);
    s.algorithm = @(x_0,lower_level_problem,upper_level_problem,sol,s,param) nonsmooth_trust_region_algorithm(lower_level_problem,upper_level_problem,sol,s,param);
    s.finalize = @(info) nonsmooth_trust_region_finalize();
end

function [sol,state,param] = nonsmooth_trust_region_initialize(x_0,lower_level_problem,upper_level_problem,param)

    state.x_n = {};
    sol = x_0;
    state.sol_history = x_0;
    state.radius = param.radius;
    state.radiusprev = param.radius;
    state.res = 1;
    state.u_history = lower_level_problem.solve(x_0,upper_level_problem.dataset);

    % Test if the min radius is defined
    param = add_default(param,'minradius',1e-4);

    % Test if using second order solver, and setting accordingly
    param = add_default(param,'use_bfgs',false);
    param = add_default(param,'use_lbfgs',false);
    param = add_default(param,'use_sr1',false);
    if param.use_bfgs==false && param.use_lbfgs==false && param.use_sr1==false
        param = add_default(param,'use_linear',true);
    end

    % Setting the hessian initialization accordingly
    if param.use_bfgs == true || param.use_lbfgs == true
        state.bfgs = 0.1*speye(size(x_0(:),1));
    elseif param.use_sr1 == true
        state.sr1 = 0.1*speye(size(x_0(:),1));
    else
      state.bfgs = sparse(size(x_0(:),1)); % Use first order model
      state.sr1 = sparse(size(x_0(:),1));
    end

end

function [sol,state] = nonsmooth_trust_region_algorithm(lower_level_problem,upper_level_problem,sol,state,param)

    % Solving the state equation (lower level solver)
    u = lower_level_problem.solve(sol,upper_level_problem.dataset);

    % Getting current cost
    cost = upper_level_problem.eval(u,sol,upper_level_problem.dataset);

    % Saving cost history
    if ~isfield(state, 'l2_cost_history')
        state.l2_cost_history = cost;
    else
        state.l2_cost_history = [state.l2_cost_history cost];
    end

    if state.radius >= param.minradius

        % Solving the Bouligand subdifferential element
        gradient_parameters.regularized_model = false;
        state.grad = upper_level_problem.gradient(u,sol,upper_level_problem.dataset,gradient_parameters);

    else

        % Solving the regularized gradient
        gradient_parameters.regularized_model = true;
        state.grad = upper_level_problem.gradient(u,sol,upper_level_problem.dataset,gradient_parameters);

    end

    % Update second order approximation matrix
    if state.radiusprev > param.minradius && state.radius < param.minradius
        if param.use_bfgs == true || param.use_lbfgs == true
            state.bfgs = 0.1*speye(size(sol(:),1));
        elseif param.use_sr1 == true
            state.sr1 = 0.1*speye(size(sol(:),1));
        end
    else
        if param.use_bfgs == true
            state = update_bfgs_approximation(sol,state);
        elseif param.use_sr1 == true
            state = update_sr1_approximation(sol,state);
        end
    end

    if param.use_bfgs == true
        step = solve_tr_subproblem(sol(:),state.grad(:),state.bfgs,state.radius);
    elseif param.use_sr1 == true
        step = solve_tr_subproblem(sol(:),state.grad(:),state.sr1,state.radius);
    else
        step = solve_tr_subproblem(sol(:),state.grad(:),0,state.radius);
    end
    step = reshape(step,size(state.grad));

    % Quality Indicator calculation
    if param.use_bfgs == true
        pred = -state.grad(:)'*step(:) - 0.5*step(:)'*state.bfgs*step(:);
    elseif param.use_sr1 == true
        pred = -state.grad(:)'*step(:) - 0.5*step(:)'*state.sr1*step(:);
    else
        pred = -state.grad(:)'*step(:);
    end
    if pred < 0 || norm(sol(:)+step(:))<1
        rho = 0;
    else
        next_u = lower_level_problem.solve(sol+step,upper_level_problem.dataset);
        next_cost = upper_level_problem.eval(next_u,sol+step,upper_level_problem.dataset);
        ared = cost-next_cost;
        rho = ared/pred;
    end

    if size(sol,1)>1 || size(sol,2)>1

        fprintf('l2_cost = %f, norm_sol = %f, norm_grad = %f, radius = %f, rho = %f, norm_step = %f',cost,norm(sol(:)),norm(state.grad(:),inf),state.radius,rho,norm(step(:)));
    else
        fprintf('l2_cost = %f, sol = %f, grad = %f, radius = %f, rho = %f, step = %f',cost,sol,state.grad,state.radius,rho,step);
    end

    if gradient_parameters.regularized_model == 1
        fprintf(' * \n');
    else
        fprintf('\n');
    end

    % Record previous radius
    state.radiusprev = state.radius;
    
    % Change size of the region
    if rho > param.eta2
        % Record previous step
        state.solprev = sol;

        % Updating solution
        sol = sol + step;
        state.radius = param.gamma2*state.radius;

        % Record previous gradient
        state.gradprev = state.grad;

        % Store information
        if size(sol,1)>1 || size(sol,2)>1
            state.sol_history = cat(3,state.sol_history,sol);
        else
            state.sol_history = [state.sol_history sol];
        end
        state.u_history = cat(3,state.u_history,u);

    elseif rho <= param.eta1
        state.radius = param.gamma1*state.radius;
    else
      %sol = sol + step;
        state.radius = param.gamma1*state.radius;
    end

end

function nonsmooth_trust_region_finalize()
    fprintf('(*) Regularized Model Evaluation\n');
end

function x = projection_linf_pos(x0,x,radius)
    % Project into the l infinity norm intersected with the positive
    % quadrant
    x_f = x0+radius;
    x_b = max(x0-radius,0);
    ind_f = find(x>x_f);
    ind_b = find(x<x_b);
    x(ind_f) = x_f(ind_f);
    x(ind_b) = x_b(ind_b);
end

function [step,directions] = solve_tr_subproblem(sol,grad,Bk,radius)
    n = length(sol);
    errtol = 1e-4;
    maxiters = 1000;
    step = sparse(n,1);
    r = -grad-Bk*step;
    rho = r'*r;
    tst = norm(r);
    terminate = errtol*norm(grad);
    it = 1;
    directions = sparse(n,1);
    hatdel = radius*(1-1.d-6);
    while((tst>terminate) && (it <= maxiters) && norm(step) <= hatdel)
        if(it==1)
            p = r;
        else
            beta = rho/rho_old;
            p = r + beta*p;
        end
        w = Bk*p;
        alpha = p'*w;
        if(alpha <0 )
            fprintf('NEGATIVE CURVATURE ');
            sol_ = sol+radius*p;
            sol_ = projection_linf_pos(sol,sol_,radius);
            step = sol_-sol;
            break;
        else
            alpha = rho/alpha;
            if norm(step+alpha*p)>radius
                sol_ = sol+alpha*p;
                sol_ = projection_linf_pos(sol,sol_,radius);
                step = sol_-sol;
                break;
            end
        end
        step = step+alpha*p;
        directions(:,it)=alpha*p;
        r = r - alpha*w;
        tst = norm(r);
        rho_old = rho;
        rho = r'*r;
        it = it+1;
    end
end

function [state] = update_bfgs_approximation(sol,state)
    if isfield(state,'gradprev')
        dsol = sol(:)-state.solprev(:);
        t = state.bfgs*dsol;
        r = state.grad(:)-state.gradprev(:);
        if norm(r) > 1e-10
            if (dsol(:)'*r(:)) < 0
                fprintf(2, '<<NEGATIVE CURVATURE - Skipping BFGS update>>')
            else
                state.bfgs=state.bfgs-1/(dsol'*t)*kron(t',t)+1/(dsol'*r)*kron(r',r);
            end
        else
            state.gradprev = state.grad;
        end
    end
end

function [state] = update_sr1_approximation(sol,state)
    if isfield(state,'gradprev')
        dsol = sol(:)-state.solprev(:);
        t = state.sr1*dsol;
        r = state.grad(:)-state.gradprev(:);
        if norm(r) > 1e-10
            u = r-t;
            if abs(dsol'*u) > 1e-10 * norm(dsol)*norm(u)
                state.sr1=state.sr1+(1/(u'*dsol))*kron(u',u);
            else
                fprintf(2, '<<SMALL DENOMINATOR - Skipping SR1 update>>')
            end
        else
            state.gradprev = state.grad;
        end
    end
end

