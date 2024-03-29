function [Xf,Xf_conf_neg,Xf_conf_pos] = tvarfima_forecast(X, N, d0,w,a,b,sigma_tv2,phi,theta,mu,sigma_arma2)
    d = calc_d(d0,w,a,b,sigma_tv2, X);
    X_fi = apply_tvfi(X, d);
    
    est = arima('Constant', mu, 'Variance', sigma_arma2, 'AR', phi, 'MA', theta);
    [Xf_arma, YMSE] = est.forecast(N, X_fi);
    
    X_fi = [X_fi; Xf_arma];
    
    for i = 1:N
        d = calc_d(d0,w,a,b,sigma_tv2, [X; 0]);
        
        L = length(X);
        
        %pi_t = pi_j(L, d(L+1));
        %Xflip = flip(X_fi);
        
        %X_fi(length(X_fi)+1) = Xf_arma(i) - sum(pi_t .* Xflip);
        
        X = apply_inv_tvfi(X_fi(1:L+1), d);
    end
    
    Xf = X(length(X)-N+1:length(X));
    
    X_fi_neg = X_fi;
    X_fi_pos = X_fi;
    
    X_fi_neg(length(X)-N+1:length(X)) = X_fi(length(X)-N+1:length(X)) - 1.92*sqrt(YMSE);
    X_fi_pos(length(X)-N+1:length(X)) = X_fi(length(X)-N+1:length(X)) + 1.92*sqrt(YMSE);
    
    X_fi_neg = apply_inv_tvfi(X_fi_neg, d);
    X_fi_pos = apply_inv_tvfi(X_fi_pos, d);
    
    Xf_conf_neg = X_fi_neg(length(X)-N+1:length(X));
    Xf_conf_pos = X_fi_pos(length(X)-N+1:length(X));
end