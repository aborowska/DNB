function [N_ess, InefFactor] = ESS(x, f)
    % f - truncation parameter
    [N,D] = size(x);
    N_ess = zeros(1,D);
    InefFactor = zeros(1,D);
    
    if (f > 0)
        for d = 1:D
            v = sum(autocorr(x(:,d),f));
            InefFactor(1,d) = (1+2*v);
            N_ess(1,d) = N/InefFactor(1,d);
        end
    else
        for d = 1:D
            [v, ~, bounds] = autocorr(x(:,d),N-1);
            L =  min(find((v < bounds(1,1)) & (v > bounds(2,1))));
            InefFactor(1,d) = (1+2*sum(v(1:L)));
            N_ess(1,d) = N/InefFactor(1,d);
        end       
    end
end 