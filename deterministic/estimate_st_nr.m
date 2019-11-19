function [Y, theta] = estimate_st_nr(X, error)

N = length(X);

wkd = [1 1 1 1 1 0 0];
wkd = repmat(wkd, 1, N-4);
wkd = [[1 1 0 0] wkd];
wkd = wkd(1:N)';

t = 1:N;

Sm = [ones(N, 1) t' sin((2*pi*(t'))./365) sin((4*pi*(t'))./365) wkd];

initial = num2cell((Sm'*Sm)\Sm'*X);
[b0,b1,c1,c3,d1] = initial{:};
initial = [b0 b1 c1 -200 + (400)*rand(1) c3 -200 + (400)*rand(1) d1];

S = @(x) x(1)*ones(N, 1) + x(2)*t' + x(3)*sin((2*pi*(t' + x(4)))./365) + x(5)*sin((4*pi*(t' + x(6)))./365) + x(7)*wkd;

L = @(x) X - S(x);
[~,R] = corrmtx(X - S(initial), length(X)-1);

Lds = @(x) L(x)'/R;

Sdtheta = @(x) [ones(N, 1) t' sin((t' + x(4))*2*pi/365) 2*x(3)*cos((t' + x(4))*2*pi/365)/365 sin((t' + x(6))*4*pi/365) 2*x(5)*cos((t' + x(6))*4*pi/365)/365 wkd];

Ldtheta = @(x) Lds(x) * Sdtheta(x);

theta = double(newton_raphson(Ldtheta, initial, error));
Y = S(theta);

end