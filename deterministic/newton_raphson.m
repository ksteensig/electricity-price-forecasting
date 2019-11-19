function x = newton_raphson(f, fj, initial, error)        
    x = initial;
    x = x - (fj(x)\f(x))';
    x_old = 0;
    
    while norm(x - x_old) > error
        x_old = x;
        x = x - (fj(x)\f(x))';
    end
    
end