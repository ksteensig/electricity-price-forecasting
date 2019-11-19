% used for sanity checking the newton_raphson procedure


f = @(x) x^2 + 2*x + 1;

x = newton_raphson(f, 0, 0.0001)