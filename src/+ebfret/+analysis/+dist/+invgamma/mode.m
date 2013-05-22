function M_theta = mode(a, b)
    % M_theta = mode(a, b)
    % 
    % Mode of an Inverse-Gamma distribution
    M_theta = b ./ (a+1);
    