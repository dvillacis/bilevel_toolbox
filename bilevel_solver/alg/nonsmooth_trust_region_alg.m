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
%   upper_level_problem: Upper Leveñ Problem struct instance
%   param: struct with the algorithm specific parameters
%     .minradius (real, default: 1e-4): model switching radius
%     .use_bfgs (boolean, default: false): use a second order model using
%                                           bfgs approximation
%     .use_lbfgs (boolean, default: false): use a second order model using
%                                           a limited memory bfgs
%     .no_bfgs (boolean, default: true): use a linear model
%
% OUTPUTS
%   sol: current state of the parameter solution
%   info: algorithm run information


function s = nonsmooth_trust_region_alg()
    s.name = 'NONSMOOTH_TRUST_REGION';
    s.initialize = @(x_0, lower_level_problem, upper_level_problem, param) nonsmooth_trust_region_initialize(x_0,lower_level_problem,upper_level_problem,param);
    s.algorithm = @(x_0,lower_level_problem,upper_level_problem,sol,s,param) nonsmooth_trust_region_algorithm(lower_level_problem,upper_level_problem,sol,s,param);
    s.finalize = @(info) nonsmooth_trust_region_finalize(info);
end

function [sol,state,param] = nonsmooth_trust_region_initialize(x_0,lower_level_problem,upper_level_problem,param)

    state.x_n = {};
    sol = x_0;
    state.sol_history = x_0;
    state.radius = param.radius;
    state.res = 1;
    state.u_history = lower_level_problem.solve(x_0,upper_level_problem.dataset);

    % Test if the min radius is defined
    param = add_default(param,'minradius',1e-4);

    % Test if using second order solver, and setting accordingly
    param = add_default(param,'use_bfgs',false);
    param = add_default(param,'use_lbfgs',false);
    param = add_default(param,'no_bfgs',true);
    
    % Setting the hessian initialization accordingly
    if param.use_bfgs == true || param.use_lbfgs == true
        state.bfgs = 0.001*speye(size(x_0(:),1));
    else
      state.bfgs = zeros(size(x_0(:),1)); % Use first order model
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

    [step,state] = tr_subproblem(sol,state,param);
    step = reshape(step,size(state.grad));
    %step = tr_generalized_cauchy(sol,state.grad,state.bfgs,state.radius,cost,param.use_bfgs);
    
    % Record previous step
    state.solprev = sol;
    state.gradprev = state.grad;
    
    % Trust Region Modification
    pred = -state.grad(:)'*step(:) - 0.5*step(:)'*state.bfgs*step(:); % TODO: Fix for limited memory bfgs
    next_u = lower_level_problem.solve(sol+step,upper_level_problem.dataset);
    next_cost = upper_level_problem.eval(next_u,sol+step,upper_level_problem.dataset);
    ared = cost-next_cost;
    rho = ared/pred;
    
    
    if size(sol,1)>1 || size(sol,2)>1
        
        fprintf('l2_cost = %f, norm_sol = %f, norm_grad = %f, radius = %f, rho = %f, norm_step = %f',cost,norm(sol),norm(state.grad(:),inf),state.radius,rho,norm(step));
    else
        fprintf('l2_cost = %f, sol = %f, grad = %f, radius = %f, rho = %f, step = %f',cost,sol,state.grad,state.radius,rho,step);
    end
    
    if gradient_parameters.regularized_model == 1
        fprintf(' * \n');
    else
        fprintf('\n');
    end
    
    % Change size of the region
    if rho > param.eta2
      sol = sol + step;
      state.radius = param.gamma2*state.radius;
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

function nonsmooth_trust_region_finalize(info)
    fprintf('(*) Regularized Model Evaluation\n');
end

function [step,state] = tr_subproblem(sol,state,param)
    % Step calculation
    %sn = -bfgs\grad; %TODO: Replace for limited memory BFGS
    %sn = get_Hg_lbgfs(grad,S,Y,hdiag);
    if isfield(state,'solprev')
        if param.use_bfgs == true
            if norm(sol-state.solprev)>1e-7
                dsol = sol(:) - state.solprev(:);
                t = state.bfgs*dsol;
                r = state.grad(:)-state.gradprev(:);
                if (dsol(:)'*r(:)) < 0
                    fprintf(2, '<<NEGATIVE CURVATURE - Skipping BFGS update>>')
                else
                    state.bfgs=state.bfgs-1/(dsol'*t)*kron(t',t)+1/(dsol'*r)*kron(r',r);
                end
            end
            sn = -state.bfgs\state.grad(:);
            predn = -state.grad(:)'*sn-0.5*sn'*state.bfgs*sn;
            if state.grad(:)'*state.bfgs*state.grad(:) <= 0 % Check curvature to see if optimizer is in the boundary or within
                t = radius/norm(state.grad(:));
            else
                t = min(norm(state.grad(:)).^2/(state.grad(:)'*state.bfgs*state.grad(:)),state.radius/(norm(state.grad(:))));
            end
            sc = -t*state.grad(:);
            predc = -state.grad(:)'*sc-0.5*sc'*state.bfgs*sc;
        elseif param.use_lbfgs == true
            %TODO implement limited memory bfgs
            error('Not yet implemented');
        else
            % Otherwise, use a linear model
            sn = -state.grad(:);
            predn = -state.grad(:)'*sn;
            t = state.radius/norm(state.grad(:));
            sc = -t*state.grad(:);
            predc = -state.grad(:)'*sc;
            
        end
    else
        % Otherwise, use a linear model
        sn = -state.grad(:);
        predn = -state.grad(:)'*sn;
        t = state.radius/norm(state.grad(:));
        sc = -t*state.grad(:);
        predc = -state.grad(:)'*sc;
    end

    % Step Selection
    if norm(sn)<=state.radius && predn >= 0.8*predc
        step = sn;
    else
        step = sc;
    end
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

function step = tr_generalized_cauchy(sol,grad,hess,radius,cost,use_bfgs)
    kubs = 0.1;
    klbs = 0.2;
    kfrd = 0.8;
    kepp = 0.25;
    tmin = 0;
    tmax = Inf;
    t = radius/norm(grad(:));
    maxit = 1000;
    it=0;
    step = zeros(size(sol));
    
    % Solve search direction
    if use_bfgs
        if find(isnan(hess))
            error('BFGS matrix became NaN');
        end
        delta=-hess\grad(:);
        delta = reshape(delta,size(grad)); % Use Newton's direction
    else
        delta=-grad; % Use gradient direction
    end

    while it < maxit
        
        sol_ = sol+t*delta;
        sol_ = projection_linf_pos(sol,sol_,radius); % Project into the positive half space intersection inf ball.
        sk = sol_-sol;
        mk = cost+grad(:)'*sk(:);
        if norm(sk(:))>radius || mk > cost + kubs*grad(:)'*sk(:)
            tmax = t;
        else
%             if norm(sk(:))>=kfrd*radius || mk >= cost + klbs*grad(:)'*sk(:)
%                 step = sk;
%                 break
%             else
%                 tmin = t;
%             end
            step = sk;
            break;
        end

        if tmax == Inf
            t = 2*t;
        else
            t = 0.5*(tmin+tmax);
        end
        it = it + 1;
    end

end


function [xi,step] = tr_subproblem_complex(grad,hess,radius)
    [m,n] = size(grad);
    B = spdiags(ones(n+1,1),1,m,n+1);
    H = B'*hess*B;
    f = [1;zeros(n,1)];
    b = zeros(m,1);
    A = [-ones(m,1),grad];
    obj = @(x) 0.5*x'*H*x + f'*x;
    nonloc = @(x) norm_constraint(x,radius);
    options = optimoptions('fmincon','Display','none');
    [x,fval] = fmincon(obj,[0;0],A,b,[],[],[],[],nonloc,options);
    xi = x(1);
    step = x(2:end);
end

function get_Hg_lbgfs(grad, S, Y, hdiag)
% This function returns the approximate inverse Hessian multiplied by the gradient, H*g
% Input
%   S:    Memory matrix (n by k) , s{i}=x{i+1}-x{i}
%   Y:    Memory matrix (n by k) , df{i}=df{i+1}-df{i}
%   g:    gradient (n by 1)
%   hdiag value of initial Hessian diagonal elements (scalar)
% Output
%   Hg    the the approximate inverse Hessian multiplied by the gradient g
% Reference
%   Nocedal, J. (1980). "Updating Quasi-Newton Matrices with Limited Storage".
%   Wiki http://en.wikipedia.org/wiki/Limited-memory_BFGS
%   two loop recursion
end
