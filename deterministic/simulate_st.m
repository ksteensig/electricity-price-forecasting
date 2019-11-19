function Y = simulate_st(params, N)
    wkd = [1 1 1 1 1 0 0];
    wkd = repmat(wkd, 1, N-4);
    wkd = [[1 1 0 0] wkd];
    wkd = wkd(1:N)';

    t = 1:N;

    S = @(x) x(1)*ones(N, 1) + x(2)*t' + x(3)*cos((2*pi*(t' - x(4)))./365) + x(5)*cos((4*pi*(t' - x(6)))./365) + x(7)*wkd;
    
    Y = S(params);
end