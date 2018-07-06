function [h, accept, A_sum] = update_h(y,h, theta, delta_h)

    T = length(y);
        
    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);
    gamma = theta(4);

    accept = 0;
    A_sum = 0;
    
    for t = 1:T    
        % Keep a record of the current h value being updated
        h_old = h(t);
        % normal  RW for h
        h(t) = h_old + delta_h*randn;
        % Calculate the log(acceptance probability):
        % Calculate the new likelihood value for the proposed move:
        % Calculate the numerator (num) and denominator (den) in turn:
        if (t == 1)
            m = mu;
            s2 = sigma2/(1-phi^2);
        else
            m = mu + phi*(h(t-1) - mu);
            s2 = sigma2;
        end
        num = -0.5*h(t) - 0.5*(y(t)^2)/exp(h(t)) + ...
              -0.5*log(s2) - 0.5*((h(t)-m)^2)/s2;
        den = -0.5*h_old - 0.5*(y(t)^2)/exp(h_old) + ...
              -0.5*log(s2) - 0.5*((h_old-m)^2)/s2;

        if (t < T)
% %                 m_2 = mu  + phi*(h(t) - mu);
% %                 num = num -0.5*log(sigma2) - 0.5*((h(t+1)-m_2)^2)/sigma2;
% %                 m_2 = mu  + phi*(h_old - mu);              
% %                 den = den -0.5*log(sigma2) - 0.5*((h(t+1)-m_2)^2)/sigma2;
%             m_2 = mu  + (1/phi)*(h(t+1) - mu);
%             num = num -0.5*log(sigma2) - 0.5*((h(t)-m_2)^2)/sigma2;
%             den = den -0.5*log(sigma2) - 0.5*((h_old-m_2)^2)/sigma2;      
            m_2 = mu  + phi*(h(t) - mu);
            m_2_old = mu  + phi*(h_old - mu);
            num = num - 0.5*(log(2*pi) + log(sigma2) + ((h(t+1)-m_2)^2)/sigma2);         
            den = den - 0.5*(log(2*pi) + log(sigma2) + ((h(t+1)-m_2_old)^2)/sigma2);          
        end
        % Proposal terms cancel since proposal distribution is symmetric.
        % All other prior terms cancel in the acceptance probability. 
        % Acceptance probability of MH step:
        A = min(1,exp(num-den));

        A_sum = A_sum + A;
        % To do the accept/reject step of the algorithm:        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
%             oldlikhood = newlikhood;     
            accept = accept+1;
        else  % Reject proposed move:
            % Na stays at current value:
            h(t) = h_old;
        end        
    end
end
