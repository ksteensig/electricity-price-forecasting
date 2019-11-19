function Xf = tvarfima_forecast(X, N, d0,w,a,b,sigma_tv2,phi,theta,mu,sigma_arma2)
    d = calc_d(d0,w,a,b,sigma_tv2, X);
    X_fi = apply_tvfi(X, d);
    
    plot(X_fi)
    
    est = arima('Constant', mu, 'Variance', sigma_arma2, 'AR', phi, 'MA', theta);
    Xf_arma = est.forecast(N, X_fi);
    
    for i = 1:N
        d = calc_d(d0,w,a,b,sigma_tv2, [X; 0]);
        
        L = length(X);
        
        pi_t = pi_j(L, d(L+1));
        Xflip = flip(X_fi);
        
        X_fi(length(X_fi)+1) = Xf_arma(i) - sum(pi_t .* Xflip);
        X = apply_inv_tvfi(X_fi, d);
    end
    
    Xf = X_fi(length(X_fi)-N:length(X_fi));
end