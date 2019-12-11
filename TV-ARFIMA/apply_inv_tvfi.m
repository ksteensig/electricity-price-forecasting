function X = apply_inv_tvfi(X_tvfi, d)
    X = X_tvfi;

    for i = 2:length(X_tvfi)
        pi_t = pi_j(i-1, d(i-1));
        Xflip = flip(X(1:i-1));
        X(i) = X_tvfi(i) - sum(pi_t .* Xflip);
    end
end