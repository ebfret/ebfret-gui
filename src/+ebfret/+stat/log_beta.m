function log_prob = log_beta(x, a, b)
    % log_prob = log_beta(x, a, b)
    %
    % Log probability distribution for a Beta distribution
    %
    % p(x | a, b) =
    %   Gamma(a+b) / (Gamma(a) * Gamma(b))
    %   x^(a-1) (1-x)^(b-1)

    x = bsxfun(@times, x, ones(size(a)));    
    a = bsxfun(@times, a, ones(size(x)));    
    b = bsxfun(@times, b, ones(size(x)));    

    log_prob = ...
        gammaln(a + b) - gammaln(a) - gammaln(b) ...
        + (a-1) .* log(x) ...
        + (b-1) .* log(1 - x);
