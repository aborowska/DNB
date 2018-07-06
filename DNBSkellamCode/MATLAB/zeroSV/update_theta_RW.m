function [theta, accept] = update_theta_RW(h, theta, hyper, delta)

    accept = zeros(1,3);
%     mu = theta(:,1);     % prior: mu ~ normpdf(c, 0, 10);
%     phi = theta(:,2);    % prior: (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5);
%     sigma2 = theta(:,3); % prior: 1/sigma2 ~ gampdf(1./s2, 5/2, 0.05/2);
%     gamma = theta(:,4);  % prior gamma ~ betapdf(gamma, 1.7, 10);

%     prior_const = [-0.5*(log(2*pi)+log(10)),  - log(beta(20, 1.5)),  2.5*log(0.025), -log(gamma(2.5))];
%     logpdf_norm = @(xx) prior_const(1,1) -0.5*(xx.^2)/10;
%     logpdf_beta = @(xx) prior_const(1,2)  + (20-1)*log(x) + (1.5-1)*log(1-xx); 
%     logpdf_invgamma = @(xx) prior_const(1,3) + prior_const(1,4) - (2.5+1)*log(xx) - 0.025./xx;
    
    % Cycle through each theta in turn 
    % and propose to update using random walk MH with uniform proposal density:
   
    for ii = 1:4
        % Keep a record of the current theta value being updated
        oldtheta = theta(ii);
        oldlikhood = loglik_h(h, theta);
        
        % Propose a new value using a RW with uniform proposal density
        theta(ii) = theta(ii) + delta(ii)*randn; 
        if ((ii == 1) || ((ii == 2) && (abs(theta(ii))<1)) || ((ii == 3) && (theta(ii)>0)))
            newlikhood = loglik_h(h, theta);

            % Calculate the log(acceptance probability):
            % Calculate the new likelihood value for the proposed move:
            % Calculate the numerator (num) and denominator (den) in turn:

            % Add in prior terms to the acceptance probability
            switch ii
                case 1
                    newprior = -0.5*(theta(ii).^2)/hyper.M;
                    oldprior = -0.5*(oldtheta.^2)/hyper.M;
                case 2
                    newprior = (hyper.P(1)-1)*log((1+theta(ii))/2) + (hyper.P(2)-1)*log((1-theta(ii))/2); 
                    oldprior = (hyper.P(1)-1)*log((1+oldtheta)/2) + (hyper.P(2)-1)*log((1-oldtheta)/2);                 
                case 3
                    newprior = - (hyper.S(1)+1)*log(theta(ii)) - hyper.S(2)./theta(ii);
                    oldprior = - (hyper.S(1)+1)*log(oldtheta) - hyper.S(2)./oldtheta;    
                otherwise  
                    newprior = (hyper.G(1)-1)*log(theta(ii)) + (hyper.G(2)-1)*log(theta(ii)); 
                    oldprior = (hyper.G(1)-1)*log(oldtheta) + (hyper.G(2)-1)*log(oldtheta);                                     
            end
            num = newlikhood + newprior;
            den = oldlikhood + oldprior;

            % All other prior terms (for other thetas) cancel in the acceptance probability.
            % Proposal terms cancel since proposal distribution is symmetric.

            % Acceptance probability of MH step:
            A = min(1,exp(num-den));
        else
            A = 0;
        end
        % To do the accept/reject step of the algorithm        
        % Accept the move with probability A:
        if (rand <= A)  % Accept the proposed move:
            % Update the log(likelihood) value:
            accept(ii) = A;
        else  % Reject proposed move:
            % theta stays at current value:
            theta(ii) = oldtheta;
        end
    end  
    
end