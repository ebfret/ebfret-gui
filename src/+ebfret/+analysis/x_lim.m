function [x_min, x_max] = x_lim(x, spread)
	% [x_min, x_max] = x_lim(x, spread)
	%
    % Returns left and right bound for observation spectrum based 
    % on a number of median deviations.
    %
    %
    % Inputs
    % ------
    % 
    % x : cell or vector
    %   Time series data
    %
    % spread : int (default: 5)
    %   Spread of means of states
    %
    % Outputs
    % -------
    %
    % x_min, x_max : scalar
    %   Limits

    if nargin < 2
        spread = 5;
    end
    if iscell(x)
        x = cat(1, x{:});
    end    
    x0 = median(x);
    left(1) = median(x(x<x0));
    right(1) = median(x(x>x0));
    for m = 2:spread
        left(m) = median(x(x<left(m-1)));
        right(m) = median(x(x>right(m-1)));
    end

    % get limits from histogram counts
    [h b] = hist(x, linspace(left(end), right(end), min(200, length(x)/10)));
    x_min = max(b(find(h > 0.01 * max(h), 1, 'first')), left(end));
    x_max = min(b(find(h > 0.01 * max(h), 1, 'last')), right(end));