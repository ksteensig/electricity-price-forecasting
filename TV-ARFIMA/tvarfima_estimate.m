function [d0,w,a,b,sigma_tv2,phi,theta,mu,sigma_arma2] = tvarfima_estimate(X, initial, lower, upper)

    XfunL = @(x) -tvfi_likelihood(x(1), x(2), x(3), x(4), x(5), X);
    dfun = @(x) calc_d(x(1), x(2), x(3), x(4), x(5), X);
    
    options = optimoptions('fmincon', 'MaxFunctionEvaluations', 10^4);
    
    tv_params = fmincon(XfunL,initial,zeros(5),zeros(5,1),zeros(5),zeros(5,1),lower,upper,[],options)

    for i = 1:16
        d = dfun(tv_params);

        U = apply_tvfi(X,d);

        %M = check_arima(U, 6, 6);

        %[row,col]=find(M==min(nonzeros(M)));

        U_arma = arima(1,0,0);
        est = U_arma.estimate(U);
        
        AR = [1 cell2mat(est.AR)];
        MA = [1 cell2mat(est.MA)];
        
        Y = filter(AR, MA, X);
        
        YfunL = @(x) -tvfi_likelihood(x(1), x(2), x(3), x(4), x(5), Y);
        
        tv_params = fmincon(YfunL,initial,zeros(5),zeros(5,1),zeros(5),zeros(5,1),lower,upper,[],options)
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