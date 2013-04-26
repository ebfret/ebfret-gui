function log_prob = log_studt(x, m, l, nu)
    % log_prob = log_studt(x, m, l, nu)
    %
    % Log probability density function for a Student t-distribution
    %
    % p(x | m, l, nu) = 
    %   Gamma((nu + 1) / 2) / Gamma(nu / 2)
    %   (l / (pi nu))^1/2 
    %   (1 + l (x - m)^2 / nu)^(-(nu + 1) / 2)
    x = bsxfun(@times, x, ones(size(m)));    
    m = bsxfun(@times, m, ones(size(x)));    
    l = bsxfun(@times, l, ones(size(x)));    
    nu = bsxfun(@times, nu, ones(size(x)));    
    log_prob = ...
        gammaln(0.5 * (nu + 1)) - gammaln(0.5 * nu) ...
        + 0.5 * log(l ./ (pi * nu)) ...
        - 0.5 * (nu + 1) .* log(1 + l .* (x - m).^2 ./ nu);
