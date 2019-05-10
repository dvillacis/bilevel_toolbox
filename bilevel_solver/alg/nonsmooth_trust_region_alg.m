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
    state.res = 1;
    state.u_history = lower_level_problem.solve(x_0,upper_level_problem.dataset);

    % Test if the min radius is defined
    param = add_default(param,'minradius',1e-4);

    % Test if using second order solver, and setting accordingly
    param = add_default(param,'use_bfgs',false);
    param = add_default(param,'use_lbfgs',false);
    param = add_default(param,'use_sr1',false);
    if param.use_bfgs==false && param.use_lbfgs==false && param.use_sr1==false
        param = add_default(param,'use_linear',false);
    end

    % Setting the hessian initialization accordingly
    if param.use_bfgs == true || param.use_lbfgs == true
        state.bfgs = 0.001*speye(size(x_0(:),1));
    elseif param.use_sr1 == true
        state.sr1 = 0.001*speye(size(x_0(:),1));
    else
      state.bfgs = zeros(size(x_0(:),1)); % Use first order model
      state.sr1 = zeros(size(x_0(:),1));
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
    if param.use_bfgs == true
        state = update_bfgs_approximation(sol,state);
    elseif param.use_sr1 == true
        state = update_sr1_approximation(sol,state);
    end

    % Solving the trust region subproblem
    if param.use_bfgs == true
        step = tr_subproblem_quadratic(sol(:),cost,state.bfgs,state.grad(:),state.radius);
    elseif param.use_sr1 == true
        step = tr_subproblem_quadratic(sol(:),cost,state.sr1,state.grad(:),state.radius);
    else
        step = tr_subproblem_linear(sol(:),cost,state.grad(:),state.radius);
    end
    step = reshape(step,size(state.grad));
    %step = tr_generalized_cauchy(sol,state.grad,state.bfgs,state.radius,cost,param.use_bfgs);

    % Quality Indicator calculation
    if param.use_bfgs == true
        pred = -state.grad(:)'*step(:) - 0.5*step(:)'*state.bfgs*step(:);
    elseif param.use_sr1 == true
        pred = -state.grad(:)'*step(:) - 0.5*step(:)'*state.sr1*step(:);
    else
        pred = -state.grad(:)'*step(:);
    end
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

function [step] = tr_subproblem_quadratic(sol,cost,Bk,grad,radius)
    %TODO: Update to a generalized cauchy point
    step = find_cauchy_point(sol,radius,grad,Bk);
%     if max(eig(Bk)) >= 0
%         search_direction = -Bk\grad;
%     else
%         search_direction = Bk\grad;
%     end
%     step = find_generalized_cauchy_point(sol,cost,search_direction,radius,grad,Bk);
    % predn = -grad'*sn-0.5*sn'*Bk*sn;
    % if grad'*Bk*grad <= 0 % Check curvature to see if optimizer is in the boundary or within
    %     t = radius/norm(grad);
    % else
    %     t = min(norm(grad).^2/(grad'*Bk*grad),radius/(norm(grad)));
    % end
    % sc = -t*grad;
    % predc = -grad'*sc-0.5*sc'*Bk*sc;
    % % Step Selection
    % if norm(sn)<=radius && predn >= 0.8*predc
    %     step = sn;
    % else
    %     step = sc;
    % end
end

function [step] = tr_subproblem_linear(sol,cost,grad,radius)
    %TODO: Update to a generalized cauchy point
    step = find_cauchy_point(sol,radius,grad);
    %search_direction = -grad;
    %step = find_generalized_cauchy_point(sol,cost,search_direction,radius,grad);
    % predn = -grad'*sn;
    % t = radius/norm(grad);
    % sc = -t*grad;
    % predc = -grad'*sc;
    % % Step Selection
    % if norm(sn)<=radius && predn >= 0.8*predc
    %     step = sn;
    % else
    %     step = sc;
    % end
end

function [step] = find_cauchy_point(sol,radius,grad,hess)
    % Check if hessian is present
    if nargin < 3
        hess = 0;
        sn = -grad;
    else
        sn1 = -hess\grad;
        sn2 = hess\grad;
    end
    
    % Dogleg strategy over a box
    predn1 = -grad'*sn1-0.5*sn1'*hess*sn1;
    predn2 = -grad'*sn2-0.5*sn2'*hess*sn2;
    if predn1 > predn2
        sn = sn1;
        predn = predn1;
    else
        sn = sn2;
        predn = predn2;
    end
    
    if grad'*hess*grad <= 0 % Check curvature to see if optimizer is in the boundary or within
        t = radius/norm(grad);
    else
        t = min(norm(grad).^2/(grad'*hess*grad),radius/(norm(grad)));
    end
    sc = -t*grad;
    predc = -grad'*sc-0.5*sc'*hess*sc;
    % Step Selection
    if norm(sn)<=radius && predn >= 0.8*predc
        step = sn;
        fprintf('SR1 ');
    else
        step = sc;
        fprintf('CAUCHY ');
    end
    
end

function [step] = find_generalized_cauchy_point(sol,cost,delta,radius,grad,hess)

    % Check if hessian is present
    if nargin < 5
        hess = 0;
    end

    t = radius/norm(grad);

    kubs = 0.2;
    tmin = 0;
    tmax = Inf;
    maxit = 1000;
    it=0;
    while it < maxit
        sol_ = sol + t*delta;
        sol_ = projection_linf_pos(sol,sol_,radius); % Project into the positive half space intersection inf ball.
        sk = sol_-sol;
        mk = cost+grad'*sk+0.5*sk'*hess*sk;
        % Checking stopping conditions
        if norm(sk)>radius || mk > cost+kubs*grad'*sk
            tmax = t;
        else
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

function [state] = update_bfgs_approximation(sol,state)
    if isfield(state,'gradprev')
        dsol = sol(:)-state.solprev(:);
        t = state.bfgs*dsol;
        r = state.grad(:)-state.gradprev(:);
        if norm(r) > 1e-8
            if (dsol(:)'*r(:)) < 0
                fprintf(2, '<<NEGATIVE CURVATURE - Skipping BFGS update>>')
            else
                state.bfgs=state.bfgs-1/(dsol'*t)*kron(t',t)+1/(dsol'*r)*kron(r',r);
            end
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
        end
    end
end


% function [step,state] = tr_subproblem(sol,state,param)
%     % Step calculation
%     %sn = -bfgs\grad; %TODO: Replace for limited memory BFGS
%     %sn = get_Hg_lbgfs(grad,S,Y,hdiag);
%     if isfield(state,'solprev')
%         if param.use_bfgs == true
%             if norm(sol-state.solprev)>1e-7
%                 dsol = sol(:) - state.solprev(:);
%                 t = state.bfgs*dsol;
%                 r = state.grad(:)-state.gradprev(:);
%                 if (dsol(:)'*r(:)) < 0
%                     fprintf(2, '<<NEGATIVE CURVATURE - Skipping BFGS update>>')
%                 else
%                     state.bfgs=state.bfgs-1/(dsol'*t)*kron(t',t)+1/(dsol'*r)*kron(r',r);
%                 end
%             end
%             sn = -state.bfgs\state.grad(:);
%             predn = -state.grad(:)'*sn-0.5*sn'*state.bfgs*sn;
%             if state.grad(:)'*state.bfgs*state.grad(:) <= 0 % Check curvature to see if optimizer is in the boundary or within
%                 t = radius/norm(state.grad(:));
%             else
%                 t = min(norm(state.grad(:)).^2/(state.grad(:)'*state.bfgs*state.grad(:)),state.radius/(norm(state.grad(:))));
%             end
%             sc = -t*state.grad(:);
%             predc = -state.grad(:)'*sc-0.5*sc'*state.bfgs*sc;
%         elseif param.use_lbfgs == true
%             %TODO implement limited memory bfgs
%             error('Not yet implemented');
%         elseif para.use_sr1 == true
%             if norm(sol-state.solprev)>1e-7
%                 dsol = sol(:) - state.solprev(:);
%                 t = state.bfgs*dsol;
%                 r = state.grad(:)-state.gradprev(:);
%                 state.sr1=state.sr1-1/(dsol'*t)*kron(t',t)+1/(dsol'*r)*kron(r',r);
%             end
%             sn = -state.sr1\state.grad(:);
%             predn = -state.grad(:)'*sn-0.5*sn'*state.sr1*sn;
%             if state.grad(:)'*state.sr1*state.grad(:) <= 0 % Check curvature to see if optimizer is in the boundary or within
%                 t = radius/norm(state.grad(:));
%             else
%                 t = min(norm(state.grad(:)).^2/(state.grad(:)'*state.sr1*state.grad(:)),state.radius/(norm(state.grad(:))));
%             end
%             sc = -t*state.grad(:);
%             predc = -state.grad(:)'*sc-0.5*sc'*state.sr1*sc;
%         else
%             % Otherwise, use a linear model
%             sn = -state.grad(:);
%             predn = -state.grad(:)'*sn;
%             t = state.radius/norm(state.grad(:));
%             sc = -t*state.grad(:);
%             predc = -state.grad(:)'*sc;
% 
%         end
%     else
%         % Otherwise, use a linear model
%         sn = -state.grad(:);
%         predn = -state.grad(:)'*sn;
%         t = state.radius/norm(state.grad(:));
%         sc = -t*state.grad(:);
%         predc = -state.grad(:)'*sc;
%     end
% 
%     % Step Selection
%     if norm(sn)<=state.radius && predn >= 0.8*predc
%         step = sn;
%     else
%         step = sc;
%     end
% end



% function step = tr_generalized_cauchy(sol,grad,hess,radius,cost,use_bfgs)
%     kubs = 0.1;
%     klbs = 0.2;
%     kfrd = 0.8;
%     kepp = 0.25;
%     tmin = 0;
%     tmax = Inf;
%     t = radius/norm(grad(:));
%     maxit = 1000;
%     it=0;
%     step = zeros(size(sol));
% 
%     % Solve search direction
%     if use_bfgs
%         if find(isnan(hess))
%             error('BFGS matrix became NaN');
%         end
%         delta=-hess\grad(:);
%         delta = reshape(delta,size(grad)); % Use Newton's direction
%     else
%         delta=-grad; % Use gradient direction
%     end
% 
%     while it < maxit
% 
%         sol_ = sol+t*delta;
%         sol_ = projection_linf_pos(sol,sol_,radius); % Project into the positive half space intersection inf ball.
%         sk = sol_-sol;
%         mk = cost+grad(:)'*sk(:);
%         if norm(sk(:))>radius || mk > cost + kubs*grad(:)'*sk(:)
%             tmax = t;
%         else
% %             if norm(sk(:))>=kfrd*radius || mk >= cost + klbs*grad(:)'*sk(:)
% %                 step = sk;
% %                 break
% %             else
% %                 tmin = t;
% %             end
%             step = sk;
%             break;
%         end
% 
%         if tmax == Inf
%             t = 2*t;
%         else
%             t = 0.5*(tmin+tmax);
%         end
%         it = it + 1;
%     end
% 
% end


% function [xi,step] = tr_subproblem_complex(grad,hess,radius)
%     [m,n] = size(grad);
%     B = spdiags(ones(n+1,1),1,m,n+1);
%     H = B'*hess*B;
%     f = [1;zeros(n,1)];
%     b = zeros(m,1);
%     A = [-ones(m,1),grad];
%     obj = @(x) 0.5*x'*H*x + f'*x;
%     nonloc = @(x) norm_constraint(x,radius);
%     options = optimoptions('fmincon','Display','none');
%     [x,fval] = fmincon(obj,[0;0],A,b,[],[],[],[],nonloc,options);
%     xi = x(1);
%     step = x(2:end);
% end
