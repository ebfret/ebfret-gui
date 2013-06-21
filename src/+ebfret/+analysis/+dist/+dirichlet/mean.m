function theta = mean(alpha)
    % theta = mean(alpha)
    % 
    % Expectation value of a Dirichlet distribution
    transposed = isrow(alpha);
    if transposed
        alpha = alpha';
    end
    theta = ebfret.normalize(alpha, 2);
    if transposed
        theta = theta';
    end
