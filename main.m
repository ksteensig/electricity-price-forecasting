%% simulate and estimate s(t)

data = csvread('elspot_prices.csv', 1, 1);

P_all = data(:,7); % DK1

N = 365*2;

P = P_all(1:N);

% estimate theta and calculate s(t) based on it
[st, thetals] = estimate_st_ls(P);

X = P - st;

% autocorr of dX indicates AR(2) and MA(6), thus ARIMA(2,1,6)
dX = filter([1,-1],1,X);

% estimate ARIMA
arima_mdl = arima(3,1,3);
em = arima_mdl.estimate(X);

st_new = simulate_st(thetals, length(X)+1);

% 1 day-ahead forecast
Xp_arma = cumsum([dX; em.forecast(1, X)]) + st_new;

%% estimate TV-ARFIMA
lower = [0.05 -inf 0.001 0.01 0];
upper = [0.495 0.08 0.15 0.1 1000];
initial = [0.25 0.1 0.05 0.05 1000];

rng default
[d0,w,a,b,sigma_tv2,phi,theta,mu,sigma_arma2] = tvarfima_estimate(X, initial, lower, upper);


%%
X_fi = arfima_estimate(X, 'FML', [0 0]);
X_arfima = X;

for i = 2:length(X)
    pi_t = pi_j(i-1, X_fi.d(1));
    Xflip = flip(X(1:i-1));
    X_arfima(i-1) = sum(pi_t .* Xflip);
end

%%
arma_mdl = arima(2,0,0);
X_arma = arma_mdl.estimate(X_arfima)
plot(X_arma.infer(X_arfima))