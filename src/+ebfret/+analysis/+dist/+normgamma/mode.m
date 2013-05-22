function [M_mu, M_lambda] = mode(m, beta, a, b)
    % [E_mu, E_lambda] = mode(m, beta, a, b)
    % 
    % Mode values for normal-gamma distributions
    M_mu = m;
    M_lambda = (a-1)./b;
    M_lambda(a<1) = nan;