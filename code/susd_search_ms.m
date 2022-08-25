function [positions,xmins,fmins,iter,s_conditions] = susd_search_ms(x0,obj_fun,z_des,max_iter,max_FEs, susd_gain, form_gain, form_dist, gain_funcs)
    % implementation of the SUSD search method
    % - x0: initial state as R^{dxN} for N agents of d dimension
    % - obj_fun: function handle for the objective function
    % - z_des: desired function minimum for stopping criteria
    % - max_iter: maximum number of iterations to search for
    % - susd_gain: gain of the susd function (defaults to 0.1)
    % - form_gain: gain of the formation control law (defaults to 0.1)
    % - form_dist: distance between agents in formation law (defaults to 5)
    % - gain_funcs: cell array of funcs to get the form_dist, susd_gain,
    % and form_gain from the minimum agent function evaluation
    if (~exist('susd_gain', 'var'))
        susd_gain = 0.1;
    end
    if (~exist('form_gain', 'var'))
        form_gain = 0.1;
    end
    if (~exist('form_dist', 'var'))
        form_dist = 5;
    end

    x = x0;
    f = @(y,x) y;

    s=0.02;
    
    dimensions = size(x0,1);
    num_agents = size(x0,2);
    swarm_count = size(x0,3);
    %swarm_num = 1;

    s_conditions = zeros(1, max_iter);
    xmins = zeros(dimensions, max_iter);
    fmins = zeros(1, max_iter);
    fmins(1) = 1e99;
    positions = zeros(dimensions,num_agents, swarm_count, max_iter);
    
    total_FEs = 0;
    
    stopped = false;
    nold = ones(dimensions,swarm_count); % temp value
    x_update=zeros(size(x));
    for iter=1:max_iter
        swarm_num = 1;
        for swarm=1:swarm_count
            % calculating minimum
            ff = obj_fun(x(:,:,swarm));
            [fmin,idx] = min(ff);
            if(exist('gain_funcs', 'var'))
               [form_dist, form_gain, susd_gain] = gain_funcs(mean(ff), x(:,idx,swarm), swarm); 
            end
            % perform a search step
            [u,z, nold(:,swarm)] = susd(obj_fun,f, x(:,:,swarm),iter,...
                                        nold(:,swarm), true);
             x_update(:,:,swarm)=susd_gain*u +form_gain*dist_form(form_dist/((iter)^(1/2)),x(:,:,swarm));      
%             x(:,:,swarm) = x(:,:,swarm)  +x_update(:,:,swarm);
%             x(:,:,swarm) = x(:,:,swarm) + susd_gain*u +...
%                            form_gain*dist_form(form_dist/iter,x(:,:,swarm));
%      
            %increment FEs (each time susd() is called, num_agents FEs occur???)
            total_FEs = total_FEs + num_agents; 
                        
            
            
            % check stopping condition
            if total_FEs >= max_FEs
                stopped = true; 
                break
            end
            
            minimum = min(z);
            if swarm == 1 || minimum < s_conditions(iter)
                s_conditions(iter) = minimum;
                if s_conditions(iter) < z_des 
                    %stopped = true; 
                    %break
                end
            end
        end
        for swarm=1:swarm_count
            % calulating minimum
            ff = obj_fun(x(:,:,swarm));
            [temp_fmin,temp_idx] = min(ff);
            if temp_fmin < fmin
                fmin = temp_fmin;
                idx = temp_idx;
                swarm_num = swarm;
            end
        end
        for swarm=1:swarm_count
            dif=(x(:,idx, swarm_num)-x(:,:,swarm));
            
            x(:,:,swarm) = x(:,:,swarm)  +(1-s)*x_update(:,:,swarm)+s*dif;
        end
        s=min(max(1/fmin*iter/max_iter*1/2,0.02), 0.4);
        
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

