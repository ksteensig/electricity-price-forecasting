function ar=check_arima(y,pp,qq)
    % pp is the maximum for p
    % qq is the maximum for q
    LOGL = zeros(pp+1,qq+1); %Initialize
    PQ = zeros(pp+1,qq+1);
    for p = 1:pp+1
        for q = 1:qq+1
            mod = arima(p-1,0,q-1);
            [fit,~,logL] = estimate(mod,y,'Display','off');
            LOGL(p,q) = logL;
            PQ(p,q) = p+q;
        end
    end
    
    LOGL = reshape(LOGL,(pp+1)*(qq+1),1);
    PQ = reshape(PQ,(pp+1)*(qq+1),1);
    [aic,~] = aicbic(LOGL,PQ+1,100);
    ar=reshape(aic,pp+1,qq+1);
    
    % the rows correspond to the AR degree (p) and the
    % columns correspond to the MA degree (q). The smallest value is best
end