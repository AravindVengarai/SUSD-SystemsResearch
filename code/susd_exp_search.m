function [positions,xmins,fmins,iter, s_conditions] = susd_exp_search(x0,obj_fun,epsilon,eps_window,max_iter,max_FEs, susd_gain, form_gain, form_dist, gain_funcs)
    % implementation of the SUSD search method
    % - x0: initial state as R^{dxN} for N agents of d dimension
    % - obj_fun: function handle for the objective function
    % - epsilon: stopping criteria epsilon (variance in moving window drops below epsilon)
    % - eps_window: moving window size for epsilon
    % - max_iter: maximum number of iterations to search for
    % - susd_gain: gain of the susd function (defaults to 0.1)
    % - form_gain: gain of the formation control law (defaults to 0)
    % - form_dist: distance between agents in formation law (defaults to 5)
    % - gain_funcs: cell array of funcs to get the form_dist, susd_gain,
    % and form_gain from the minimum agent function evaluation
    if (~exist('susd_gain', 'var'))
        susd_gain = 5;
    end
    if (~exist('form_gain', 'var'))
        form_gain = 0;
    end
    if (~exist('form_dist', 'var'))
        form_dist = 5;
    end

    x = x0;
    f = @(y,x) 1-exp(min(y)-y);

    dimensions = size(x0,1);
    num_agents = size(x0,2);
    swarm_count = size(x0,3);

    s_conditions = zeros(1, max_iter);
    xmins = zeros(dimensions, max_iter);
    fmins = zeros(1, max_iter);
    fmins(1) = 1e99;
    positions = zeros(dimensions,num_agents, swarm_count, max_iter);
    
    total_FEs = 0;
    
    cond = zeros(1,eps_window, swarm_count);
    stopped = false;
    nold = zeros(dimensions,swarm_count); % temp value
    for iter=1:max_iter
        swarm_num = 1;
        for swarm=1:swarm_count
            % calculating minimum
            ff = obj_fun(x(:,:,swarm));
            [fmin,idx] = min(ff);
            if(exist('gain_funcs', 'var'))
               [form_dist, form_gain, susd_gain] = gain_funcs(mean(ff), x(:,idx,swarm), swarm); 
            end
            range_without_min = [1:idx-1 idx+1:size(x,2)];
            % perform a search step
            [u,z, nold(:,swarm)] = susd(obj_fun,f, x(:,:,swarm),iter,...
                                        nold(:,swarm), false);
            x(:,:,swarm) =  x(:,:,swarm) + susd_gain*u;
            %apply formation to all agents except min
            form_change = dist_form(form_dist, x(:,:,swarm));
            x(:,range_without_min, swarm) =  x(:,range_without_min, swarm)+...
                form_gain*form_change(range_without_min);
            %increment FEs (each time susd() is called, num_agents FEs occur???)
            total_FEs = total_FEs + num_agents; 
                         
            % calulating minimum
            ff = obj_fun(x(:,:,swarm));
            [temp_fmin,temp_idx] = min(ff);
            if temp_fmin < fmin
                fmin = temp_fmin;
                idx = temp_idx;
                swarm_num = swarm;
            end
                         
            % check stopping condition
            if total_FEs >= max_FEs
                stopped = true; 
                break
            end
            
            cond(:,:,swarm) = circshift(cond(:,:,swarm),-1);
            cond(1,1,swarm) = min(z);
            s_condition = var(cond(:,:,swarm));
            if swarm == 1 || s_condition < s_conditions(iter)
                s_conditions(iter) = s_condition;
                if iter > eps_window && s_conditions(iter) < epsilon 
                    stopped = true; 
                end
            end
        end
        xmins(:,iter) = x(:,idx, swarm_num);
        fmins(iter) = min(fmin,min(setdiff(fmins,min(fmins))));
        positions(:,:,:,iter) = x;
        if stopped; break; end
    end
    
    xmins = xmins(:,1:iter);
    fmins = fmins(1:iter);

    s_conditions = s_conditions(1:iter);
    positions = positions(:,:,:,1:iter);
end