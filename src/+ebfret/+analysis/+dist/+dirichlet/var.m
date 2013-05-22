function var_theta = var(alpha)
    % var_theta = var(alpha)
    % 
    % Variance of a dirichlet distribution
    transposed = isrow(alpha);
    if transposed
        alpha = alpha';
    end
    Alpha = bsxfun(@times, sum(alpha, 2), ones(1,size(alpha,2)));
    var_theta = alpha .* (Alpha - alpha) ./ (Alpha.^2 .* (Alpha + 1));
    if transposed
        var_theta = var_theta';
    end
