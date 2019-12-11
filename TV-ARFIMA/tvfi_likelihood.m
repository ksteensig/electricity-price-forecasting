function [L] = tvfi_likelihood(d0, w, a, b, sigma_tv2, X)
        d = calc_d(d0, w, a, b, sigma_tv2, X);
        
        %L = zeros(length(X)-1,1);

        L = -0.5*log(sigma_tv2) - 1/(2*sigma_tv2) * X(1)^2;
        
        for t = 2:length(X)-1
            pi_t = pi_j(t-1, d(t));
            Xflip = flip(X(1:t-1));
            
            L = L - 0.5*log(sigma_tv2) - 1/(2*sigma_tv2) * (X(t) + sum(pi_t .* Xflip))^2;
            %L(t) = (X(t) + sum(pi_t .* Xflip));
        end
end