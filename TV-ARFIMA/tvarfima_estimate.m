function [d0,w,a,b,sigma_tv2,phi,theta,mu,sigma_arma2] = tvarfima_estimate(X, initial, lower, upper)
    fun = @(x) -tvfi_likelihood(x(1), x(2), x(3), x(4), x(5), X)/length(X);
    dfun = @(x) calc_d(x(1), x(2), x(3), x(4), x(5), X);
    
    tv_params = [];
    %options = optimoptions('fmincon','Display','off');

    for i = 1:1
        tv_params = fmincon(fun,initial,zeros(5),zeros(5,1),zeros(5),zeros(5,1),lower,upper);

        d = dfun(tv_params);

        X_tvfi = apply_tvfi(X,d);

        M = check_arima(X_tvfi, 6, 4);

        [row,col]=find(M==min(nonzeros(M)));

        X_arma = arima(row,0,col);
        est = X_arma.estimate(X_tvfi);
        
        mu = est.Constant;
        sigma_arma2 = est.Variance;
        phi = est.AR;
        theta = est.MA;
    end
    
    d0 = tv_params(1);
    w = tv_params(2);
    a = tv_params(3);
    b = tv_params(4);
    sigma_tv2 = tv_params(5);
end