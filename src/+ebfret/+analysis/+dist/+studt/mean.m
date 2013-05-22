function V_theta = mean(l, nu)
    % V_theta = var(l, nu)
    % 
    % Variance of a Student T
    V_theta = nu ./ ((nu-2) .* l);
    V_theta(nu<=2) = nan;
