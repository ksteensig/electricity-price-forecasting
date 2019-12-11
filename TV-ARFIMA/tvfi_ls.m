function [L] = tvfi_ls(d0, w, a, b, sigma_tv2, X)
        d = calc_d(d0, w, a, b, sigma_tv2, X);
        
        L = zeros(length(X)-1,1);

        L(1) = X(1);
        
        for t = 2:length(X)-1
            pi_t = pi_j(t-1, d(t));
            Xflip = flip(X(1:t-1));
            
            L(t) = (X(t) + sum(pi_t .* Xflip));
        end
end