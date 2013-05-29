function bins = get_bins(x, num_bins, threshold)
    % bins = get_bins(x, num_bins, threshold)
    % 
    % Estimates a set of bin positions that is less sensitive to
    % outliers than Matlab's histogram function
    x = x(:);
    % ignore inf and nan points
    x = x(isfinite(x));
    if nargin < 2
        num_bins = max(10, min(200, round(length(x)/10))); 
    end
    if nargin < 03
        threshold = 0.001;
    end
    % average number of observations per bin
    dens = length(x) / num_bins;
    % calculate median deviations from median until 
    % fraction of outliers drops below threshold
    xm = median(x);
    l = find(x<xm);
    while (length(l) / length(x)) > threshold
        xl = median(x(l));
        l = find(x < xl);
    end
    r = find(x>xm);
    while (length(r) / length(x)) > threshold
        xr = median(x(r));
        r = find(x > xr);
    end
    bins = linspace(xl, xr, num_bins);