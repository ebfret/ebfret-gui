function [E_mu, E_lambda] = mean(m, beta, a, b)
    % [E_mu, E_lambda] = mean(m, beta, a, b)
    % 
    % Expectation values for normal-gamma distributions
    E_mu = m;
    E_lambda = a./b;
