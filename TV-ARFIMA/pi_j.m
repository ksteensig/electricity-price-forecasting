function [pi_t] = pi_j(j, dt)
    pi_t = ones(j,1);
    
    pi_t(1) = -dt;
    
    for i = 2:j
        pi_t(i) = pi_t(i-1) * (i - 1 - dt)/i;
    end
    
    %lower = (1:j).^-dt;
    %upper = (0:(j-1)).^-dt;
    %average = (lower+upper)/2;
    
    %pi_t = lower ./ (gamma(-dt) .* (1:j));
    %pi_t = pi_t';
end