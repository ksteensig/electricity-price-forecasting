d0 = 0.25;
w = 0.05;
a = 0.1;
b = 0.1;
sigma_tv2 = 1000;

N = 200;

lower = [0.001 -0.2 -1 -0.999 0];
upper = [0.499 0.2 1 0.999 inf];
initial = [d0 w a b sigma_tv2];

arma = arima('Constant', 3, 'AR', {-0.5}, 'MA', {0.5}, 'Variance', 1000);
U = arma.simulate(N);
X = U;

for i=2:N
    d = calc_d(d0, w, a, b, sigma_tv2, X(1:i-1));
    pi_t = pi_j(i-1, d(i-1));
    Xflip = flip(X(1:i-1));
    X(i) = U(i) - sum(pi_t .* Xflip);
end

[d0,w,a,b,sigma_tv2,phi,theta,mu,sigma_arma2] = tvarfima_estimate(X, initial, lower, upper);
d = calc_d(d0,w,a,b,sigma_tv2, X);
rmse = sqrt(sum((apply_tvfi(X,d)-U).^2)/N)

%%
index = 1;
rmse = zeros(19,1);

for i = [10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000]
    [d0,w,a,b,sigma_tv2,phi,theta,mu,sigma_arma2] = tvarfima_estimate(X, initial, lower, upper);
    d = calc_d(d0,w,a,b,sigma_tv2, X(1:i));
    rmse(index) = sqrt(sum((apply_tvfi(X(1:i),d)-U(1:i)).^2)/N);
    index = index+1;
end

%%
d = calc_d(d0,w,a,b,sigma_tv2, X);
plot(apply_tvfi(X,d))
hold on
plot(U)
rmse = sqrt(sum((apply_tvfi(X,d)-U).^2)/N)


%%
d = calc_d(0.3404, 0, 0, 1, 0, X);
X = apply_tvfi(P()-st,ones(729,1)*0.3404);

arfima = arima('Constant', 0, 'AR', {0.2545213}, 'Variance', 1620.023)

%arfima = estimate(arima(1,1,1), X);

em = [P-st;apply_inv_tvfi(arfima.forecast(300, X), ones(1030,1)*0.3404)];