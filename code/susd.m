function [u,z, nold] = susd(obj_fun,f,x,iter, nold, normalize)
    if (~exist('normalize', 'var'))
        normalize = false;
    end

    % implementation of SUSD with function mapping f
    
    % compute the eigenvectors/SUSD direction
    umean = mean(x');
    R_u = x' - repmat(umean, size(x,2), 1);
    cov = R_u'*R_u;
    [eig_vec, ~] = eig(cov);
    n = eig_vec(:,1);
    
    
    % get function evaluation
    z = obj_fun(x);
    z_ = f(z,x);
    
    if(normalize)
        z_ = z_/(max(z_)+1e-10).*(1-exp(-z_));
        %z_ = 1-exp(min(z_)-z_); 
    end
    
    
    % check for switching
    %if isempty(nold); nold = n; end
    if iter == 1
        nold = n;
        az = z_*R_u*n;
        if az>0
            sw = -1;
        else
            sw = 1;
        end
        n = sw*n;
    end
    if sign(n'*nold) ~= 0
        n = n*sign(n'*nold);
    end
    nold = n;

%     sum_val = sum(z'.*R_u*n);
%     if sum_val >= 0 
%        n = -n; 
%     end
    % compute input
    u = z_.*n;
end

