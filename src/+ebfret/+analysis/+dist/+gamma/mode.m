function theta = mode(a, b)
    % theta = mode(a, b)
    % 
    % Mode of a Gamma distribution
    theta = (a-1) ./ b;
    theta(theta<1) = nan;