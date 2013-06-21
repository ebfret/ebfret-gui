function lines = scale(lines, s)
    % lines = scale(lines, s)
    % 
    % Rescales a set of lines by multiplication with s. 
    if isscalar(s)
        s = repmat(s, size(lines));
    end
    for l = 1:length(lines(:))
        lines(l).ydata = bsxfun(@times, s(l,:), lines(l).ydata);
    end
