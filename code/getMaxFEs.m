function maxFEs = getMaxFEs(dimension)
    if dimension == 2
        maxFEs = 10000;
    elseif dimension == 10
        maxFEs = 200000;
%     elseif dimension == 15
%         maxFEs = 3000000;
    elseif dimension == 20
        maxFEs = 1000000;
    else
        error("invalid dimension for CEC2021")
    end