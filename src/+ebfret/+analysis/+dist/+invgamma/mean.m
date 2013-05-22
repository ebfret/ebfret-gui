function theta = mean(a, b)
    % theta = mean(a, b)
    % 
    % Expectation value of a Inverse-Gamma distribution
    theta = b ./ (a-1);
    theta(a<=1) = nan;