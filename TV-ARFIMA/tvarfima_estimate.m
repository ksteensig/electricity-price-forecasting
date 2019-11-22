function [d0,w,a,b,sigma_tv2,phi,theta,mu,sigma_arma2] = tvarfima_estimate(X, initial, lower, upper)

    Xfun = @(x) -tvfi_likelihood(x(1), x(2), x(3), x(4), x(5), X)/length(X);

    dfun = @(x) calc_d(x(1), x(2), x(3), x(4), x(5), X);
    
    %options = optimoptions('fmincon','Display','off');
    
    tv_params = fmincon(Xfun,initial,zeros(5),zeros(5,1),zeros(5),zeros(5,1),lower,upper);

    for i = 1:4
        d = dfun(tv_params);

        U = apply_tvfi(X,d);

        %M = check_arima(U, 5, 5);

        %[row,col]=find(M==min(nonzeros(M)));

        U_arma = arima(1,0,0);
        %est = U_arma.estimate(U, 'Display', 'off');
        est = U_arma.estimate(U);
        
        AR = cell2mat(est.AR);
        MA = cell2mat(est.MA);
        
        if isempty(est.AR)
            AR = 1;
        end
        
        if isempty(est.MA)
            MA = 1;
        end
        
        Y = filter(AR, MA, X);
        
        Yfun = @(x) -tvfi_likelihood(x(1), x(2), x(3), x(4), x(5), Y)/length(Y);
        
        tv_params = fmincon(Yfun,initial,zeros(5),zeros(5,1),zeros(5),zeros(5,1),lower,upper);
        
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