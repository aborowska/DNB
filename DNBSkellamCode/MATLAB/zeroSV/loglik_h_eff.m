function loglik = loglik_h_eff(h, theta)

    T = length(h);
    
    mu = theta(:,1);
    phi = theta(:,2);
    sigma2 = theta(:,3);
    
%     a1 = mu;
%     P1 = sigma2./(1-phi.^2);
    
    loglik = zeros(1,T);
        
    loglik(1,1) =  ((h(1) - mu).^2).*(1-phi.^2);
    for t = 2:T       
        loglik(1,t) = (h(t) - mu - phi*(h(t-1) - mu))^2;    
    end
    loglik = loglik/sigma2;
    loglik(1,1) =  loglik(1,1) - log(1-phi^2);
    
    loglik = loglik  + log(sigma2) + log(2*pi);
    loglik = -0.5*loglik;    
end