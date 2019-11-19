function [pi_t] = pi_j(j, dt)
    pi_t = ones(j,1);
    
    pi_t(1) = 1 -1 - dt;
    
    for i = 2:j
        pi_t(i) = pi_t(i-1) * (i - 1 - dt)/i;
    end
end