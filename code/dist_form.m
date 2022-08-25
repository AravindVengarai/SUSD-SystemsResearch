function [u] = dist_form(d0,x,normalize_distance)
    if (~exist('normalize_distance', 'var'))
        normalize_distance = false;
    end
    % implementation of distance based formation controller
    u = zeros(size(x,1),size(x,2));
    num_agents = size(x,2);
    
    % compute the formation term
    for i=1:num_agents
        x_ = x; 
        x_(:,i) = [];
        dif = x_ - x(:,i);
        d = vecnorm(dif);
        if normalize_distance
            u(:,i) = sum((d-d0).*dif./(d.^2),2);
        else
            u(:,i) = sum((d-d0).*dif,2);
        end
    end
end

