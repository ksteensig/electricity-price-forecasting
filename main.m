%% simulate and estimate s(t)

data = csvread('elspot_prices.csv', 1, 1);

P_all = data(:,7); % DK1

N = 365*2;

%% estimate ARIMA
arima_mdl = arima(2,1,2);
arima_fx = zeros(length(P_all)-N,2);

P = P_all(1:N);
[st, thetals] = estimate_st_ls(P);
X = P - st;

em = arima_mdl.estimate(X);

for i = N+1:length(P_all)
    P = P_all(1:i-1);    
    st = simulate_st(thetals, i);
    X = P - st(1:i-1);
    
    [fx, ymse] = em.forecast(1, X);
    
    arima_fx(i-N,1:2) = [fx ymse];
end

fxn = arima_fx(:,1) - 1.96*sqrt(arima_fx(:,2)) + st(end);
fxp = arima_fx(:,1) + 1.96*sqrt(arima_fx(:,2)) + st(end);
fx = arima_fx(:,1) + st(end);

rmse_arima = sqrt(sum((P_all(N+1:length(P_all)) - fx).^2)/(length(P_all)-N));

%% estimate TV-ARFIMA
lower = [0.05 -inf 0.001 0.01 0];
upper = [0.495 0.08 0.15 0.1 1000];
initial = [0.25 0.1 0.05 0.05 1000];

tvarfima_fx = zeros(length(P_all)-N,3);

X = P_all(1:730);

[d0,w,a,b,sigma_tv2,phi,theta,mu,sigma_arma2] = tvarfima_estimate(X, initial, lower, upper);
[~, thetals] = estimate_st_ls(P);

for i = N+1:length(P_all)
    P = P_all(1:i-1);
    
    st = simulate_st(thetals, i);
    X = P - st(1:i-1);
    
    [fx, fxn, fxp] = tvarfima_forecast(X, 1, d0,w,a,b,sigma_tv2,phi,theta,mu,sigma_arma2);
    
    tvarfima_fx(i-N,1:3) = [fx fxn fxp] + st(end);
end

rmse_tvarfima = sqrt(sum(((P_all(N+1:length(P_all)) - tvarfima_fx(:,1)).^2)/(length(P_all)-N)));

%%
X_fi = arfima_estimate(X, 'FML', [0 0]);

%% estimate ARFIMA
arima_mdl = arima(1,0,1);
arima_fx = zeros(length(P_all)-N,3);

for i = N+1:length(P_all)
    P = P_all(1:i-1);
    
    st = simulate_st(thetals, i);
    X = P - st(1:i-1);
    
    [fx, fxn, fxp] = tvarfima_forecast(X, 1, d0,w,a,b,sigma_tv2,phi,theta,mu,sigma_arma2);
    
    tvarfima_fx(i-N,1:3) = [fx fxn fxp] + st(end);
end

fxn = arima_fx(:,1) - 1.96*sqrt(arima_fx(:,2)) + arima_fx(:,3);
fxp = arima_fx(:,1) + 1.96*sqrt(arima_fx(:,2)) + arima_fx(:,3);
fx = arima_fx(:,1) + arima_fx(:,3);

rmse_arima = sqrt(sum((P_all(N+1:length(P_all)) - fx).^2)/(length(P_all)-N));

%%
arma_mdl = arima(1,0,1);
X_arma = arma_mdl.estimate(X_arfima)
plot(X_arma.infer(X_arfima))