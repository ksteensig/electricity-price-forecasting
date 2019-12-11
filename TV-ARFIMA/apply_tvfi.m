function X_tvfi = apply_tvfi(X, d)
    X_tvfi = X;

    for i = 2:length(X)
        pi_t = pi_j(i-1, d(i-1));
        Xflip = flip(X(1:i-1));
        X_tvfi(i) = X(i) + sum(pi_t .* Xflip);
    end
end