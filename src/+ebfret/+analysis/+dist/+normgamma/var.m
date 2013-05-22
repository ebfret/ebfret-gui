function [V_mu, V_lambda] = var(m, beta, a, b)
    % [E_mu, E_lambda] = mode(m, beta, a, b)
    % 
    % Mode values for normal-gamma distributions
    V_mu = b ./ (beta .* (a-1));
    V_lambda = a ./ b.^2;
