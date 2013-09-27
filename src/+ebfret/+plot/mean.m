function lines = mean(lines, dim)
    % lines = mean(lines, dim)
    % 
    % Take mean of ydata along specified dimension
    for l = 1:length(lines(:))
        lines(l).ydata = mean(lines(l).ydata, 2);
    end
