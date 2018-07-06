function [ll, rec] = loglik_h(h, theta)

    T = length(h);
%     M = size(theta,1);
    
    mu = theta(:,1);
    phi = theta(:,2);
    sigma2 = theta(:,3);
     
        

    rec =  ((h(1) - mu)^2)*(1-phi.^2) ;
    for t = 2:T
        rec = rec + (h(t) - mu - phi*(h(t-1) - mu))^2;    
    end
    ll  = rec/sigma2  - log(1-phi^2) + T*log(sigma2) + T*log(2*pi);
    ll  = -0.5*ll;
     
end