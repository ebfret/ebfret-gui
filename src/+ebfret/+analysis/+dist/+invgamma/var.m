function V_theta = var(a, b)
    % V_theta = var(a, b)
    % 
    % Variance of a Inverse-Gamma distribution
    V_theta = b.^2 ./ ((a-1) .* (a-2));
    V_theta(a<=2) = nan;