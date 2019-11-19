data = csvread('elspot_prices.csv', 1, 1);

X = data(:,7);
[Yls, thetals] = estimate_st_ls(X);

Z = X-Yls;
%dZ = filter([-1,1],1,Z);

em = estimate(arima(5,1,4), Z(1:2*365));
e = infer(em,Z(1:2*365));

dfun = @(x) calc_d(x(1), x(2), x(3), x(4), x(5), e');
d = dfun([0.2725   -0.04   0.1105    0.0839  999.9987]);

Y = zeros(2*365, 1);

for i = 2:2*365
    pi_t = pi_j(i-1, d(i-1));
    Xflip = flip(e(1:i-1))';
    Y(i-1) = sum(pi_t .* Xflip);
end

%plot(Z(1:2*365))
plot(e)
hold on
%plot(Y)
plot(e(1:2*365)+Y)


% together they imply stationarity
adftest(e(1:2*365)+Y)
kpsstest(e(1:2*365)+Y)