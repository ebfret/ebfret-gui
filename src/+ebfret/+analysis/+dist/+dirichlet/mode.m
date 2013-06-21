function theta = mode(alpha)
    % theta = mode(alpha)
    % 
    % Mode of a dirichlet distribution
    transposed = isrow(alpha);
    if transposed
        alpha = alpha';
    end
    theta = ebfret.normalize(alpha - 1, 2);
    k = find(sum(alpha < 1, 2));
    theta(k, :) = nan;
    if transposed
        theta = theta';
    end
