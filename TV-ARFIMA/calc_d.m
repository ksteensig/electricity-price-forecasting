function [d] = calc_d(d0, w, a, b, v, X)
    g = zeros(1,length(X));
    s = zeros(1,length(X));
    
    nu_j = @(j, d) pi_j(length(j), d).*(-psi(j-d) + psi(1-d) + 1/d);
    
    % derivative of d_to_g with respect to d
    hdot = @(d) -((2*d - 1)*((4*d)/(2*d - 1)^2 - 2/(2*d - 1)))/(2*d);
    
    d_to_g = @(d) log(-2*d./(2*d - 1));
    
    % 0.4999 is used instead of 0.5 to ensure the result is always valid
    g_to_d = @(g) 0.01 + (0.499 - 0.01)./(1 + exp(-g));
    
    g0 = d_to_g(d0);
    g(1) = w + b*g0; % s(1) = 0 by definition
    
    for t = 2:length(X)-1
        dt = g_to_d(g(t));
        pi_t = pi_j(t-1, dt);
        nu_t = nu_j((1:t-1)', dt);
        Xflip = flip(X(1:t-1));

        %St = v*sum(nu_t)^(-2);
        Delt = -(X(t) + sum(pi_t .* Xflip))*sum(nu_t .* Xflip)/v;
        
        hdot_t = hdot(dt);
        
        %s(t) = St * hdot_t * Delt;
        s(t) = hdot_t * Delt;
        
        % explanation of this in the GAS R package docs
        g(t+1) = w + a*s(t) + b*g(t);
    end
    
    d = g_to_d(g)';
end