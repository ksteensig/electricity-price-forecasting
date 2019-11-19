function [Y, theta] = estimate_st_ls(X)

N = length(X);

wkd = [1 1 1 1 1 0 0];
wkd = repmat(wkd, 1, N-4);
wkd = [[1 1 0 0] wkd];
wkd = wkd(1:N)';

t = 1:N;

Sm = [ones(N,1) t' cos((2*pi*(t'))./365) sin((2*pi*(t'))./365) cos((4*pi*(t'))./365) sin((4*pi*(t'))./365) wkd];
S = @(x) x(1)*ones(N, 1) + x(2)*t' + x(3)*cos((2*pi*(t' - x(4)))./365) + x(5)*cos((4*pi*(t' - x(6)))./365) + x(7)*wkd;

initial = num2cell((Sm'*Sm)\Sm'*X);
[b0,b1,c1,c2,c3,c4,d1] = initial{:};
initial = [b0 b1 c1 c2 c3 c4 d1];

[~,R] = corrmtx(X - S(initial), length(X)-1);

f = (Sm'/R*Sm)\Sm'/R*X;

theta = [ f(1); f(2); sqrt(f(3)^2 + f(4)^2); 365*atan(f(4)/f(3))/(2*pi); sqrt(f(5)^2 + f(6)^2); 365*atan(f(6)/f(5))/(2*pi); f(7)];
Y = S(theta);

end