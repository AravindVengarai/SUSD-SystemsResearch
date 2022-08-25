function positions = generateStartingPositions(dimension, min, max, num_agents, seed, num_swarms, use_form_law, distance)
    if (~exist('num_swarms', 'var'))
        num_swarms = 1;
    end    
    if (~exist('use_form_law', 'var'))
        use_form_law = false;
    end
    if (~exist('distance', 'var'))
        distance = 0.1;
    end    
    positions = zeros(dimension,num_agents, num_swarms);
    rng(seed,'twister');
        
    for i=1:num_swarms
        for j = 1:dimension
            pos = unifrnd(min, max);
            for k=1:num_agents
                if pos > 0
                    %min should be less than 0, so this should put all
                    %agents slightly lower than the pos randomly;
                    x0 = pos + unifrnd(min, 0)*.01;
                else
                    x0 = pos + unifrnd(0, max)*.01; 
                end
                positions(j,k,i) = x0;
            end
        end
        if use_form_law
            for a = 1:5
                    positions(:,:,i) = positions(:,:,i) + 0.1*dist_form(distance,positions(:,:,i), true);
            end
        end
    end
end