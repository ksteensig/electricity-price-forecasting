%% simulate and estimate s(t)

data = csvread('elspot_prices.csv', 1, 1);

P_all = data(:,7); % DK1

N = 2*365;

P_all = P_all(1:2*N);

flength = length(P_all) - N; % forecast length

fplot = 731:761;

%% estimate ARIMA
arima_mdl = arima(1,1,1);
arima_fx = zeros(length(P_all)-N,3);

P = P_all(1:N);
[st, thetals] = estimate_st_ls(P);
X = P - st;

em = arima_mdl.estimate(X);

for i = N+1:length(P_all)
    P = P_all(1:i-1);
    st = simulate_st(thetals, i);
    X = P - st(1:i-1);
    
    [fx, ymse] = em.forecast(1, X);
    
    fxn = fx - 1.96*sqrt(ymse);
    fxp = fx + 1.96*sqrt(ymse);
    
    arima_fx(i-N,1:3) = [fx fxn fxp] + st(end);
end

rmse_arima = sqrt(sum((P_all(N+1:end) - arima_fx(:,1)).^2)/flength);

figure('Name', 'ARIMA')
plot(731:761, P_all(731:761))
hold on
plot(731:761, arima_fx(1:31,1))
plot(731:761,arima_fx(1:31,2),'r:','LineWidth',2)
plot(731:761,arima_fx(1:31,3),'r:','LineWidth',2)
xlabel('Day')
ylabel('DKK/MWh')
legend('Real price', 'Forecasted price', '95% confidence')

%% estimate TV-ARFIMA
%lower = [0.05 -0.1 -0.2 -0.999 1];
%upper = [0.495 0.1 0.2 0.999 inf];
initial = [0.34 0 0.005 0.05 1];

lower = [0.05 -inf -inf -0.999 1];
upper = [0.495 inf inf 0.999 inf];
%initial = [0.25 0 0 0 1000];

%lower = [0.0001 -inf 0 -0.999 1];
%upper = [0.4999 inf inf 0.999 inf];
%initial = [0.4 0 0 1 1000];

tvarfima_fx = zeros(length(P_all)-N,3);

P = P_all(1:N);
[st, thetals] = estimate_st_ls(P);
X = P - st;

[d0,w,a,b,sigma_tv2,phi,theta,mu,sigma_arma2] = tvarfima_estimate(X, initial, lower, upper);

for i = N+1:length(P_all)
    P = P_all(1:i-1);
    
    st = simulate_st(thetals, i);
    X = P - st(1:i-1);
    
    [fx, fxn, fxp] = tvarfima_forecast(X, 1, d0,w,a,b,sigma_tv2,phi,theta,mu,sigma_arma2);
    
    tvarfima_fx(i-N,1:3) = [fx fxn fxp] + st(end);
end

rmse_tvarfima = sqrt(sum((P_all(N+1:end) - tvarfima_fx(:,1)).^2)/flength)

%%
figure('Name', 'TV-ARFIMA')
plot(731:761, P_all(731:761))
hold on
plot(731:761, tvarfima_fx(1:31,1))
plot(731:761,tvarfima_fx(1:31,2),'r:','LineWidth',2)
plot(731:761,tvarfima_fx(1:31,3),'r:','LineWidth',2)
xlabel('Day')
ylabel('DKK/MWh')
legend('Real price', 'Forecasted price', '95% confidence')


%% estimate ARFIMA
arfima_fx = zeros(length(P_all)-N,1);

P = P_all(1:N);
[st, thetals] = estimate_st_ls(P);
X = P - st;

em = arfima_estimate(X, 'FWHI', [6 5]);

for i = N+1:length(P_all)
    P = P_all(1:i-1);
    
    st = simulate_st(thetals, i);
    X = P - st(1:i-1);
    
    fx = arfima_forecast(X, 1, em.d(1), em.AR, em.MA, em.mean, em.sigma2);
    
    arfima_fx(i-N) = fx + st(end);
end

%rmse_arfima = sqrt(sum(((P_all(N+1:length(P_all)) - arfima_fx).^2)/(length(P_all)-N)));

plot(731:761, P_all(731:761))
hold on
plot(731:761, arfima_fx(1:31,1))
%plot(731:761,arfima_fx(1:31,2),'r:','LineWidth',2)
%plot(731:761,arfima_fx(1:31,3),'r:','LineWidth',2)
xlabel('Day')
ylabel('DKK/MWh')
legend('Real price', 'Forecasted price')