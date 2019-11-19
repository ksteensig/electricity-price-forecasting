function [L] = tvfi_likelihood(d0, w, a, b, v, X)
        d = calc_d(d0, w, a, b, v, X);

        L = -0.5 * log(v) - 1/(2*v) * X(1)^2;
        
        for t = 2:length(X)-1
            pi_t = pi_j(t-1, d(t));
            Xflip = flip(X(1:t-1));
            
            L = L -0.5*log(v) - 1/(2*v) * (X(t) + sum(pi_t .* Xflip))^2;
            d(t);
        end
end