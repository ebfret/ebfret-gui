function log_prob = log_gamma(x, a, b)
    % log_prob = log_gamma(x, a, b)
    %
    % Log probability distribution for a Gamma distribution
    %
    % p(x | a, b) =
    %   b^a / Gamma(a) x^(a-1) exp(-b x)

    x = bsxfun(@times, x, ones(size(a)));    
    a = bsxfun(@times, a, ones(size(x)));    
    b = bsxfun(@times, b, ones(size(x)));    

    log_prob = ...
        a .* log(b) ...
        - gammaln(a) ...
        + (a - 1) .* log(x) ...
        - b .* x;
