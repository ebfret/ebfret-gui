function u = guess_prior(x, K, spread, strength)
	% u = guess_prior(x, K, spread)
	%
    % Returns initial guess for prior based on observation histogram.
    %
    %
    % Inputs
    % ------
    % 
    % x : cell
    %   Time series data
    %
    % spread : int (default: 5)
    %   Spread of means of states
    %
    % strength : scalar (default: 1)
    %   Scale prior counts by this factor
    %
    %
    % Outputs
    % -------
    %
    % u : struct
    %   Prior parameters 
    if nargin < 3
        spread = 5;
    end
    if nargin < 4
        strength = 1;
    end
    T = cellfun(@length, x);
    x = cat(1, x{:});
    [x_min, x_max] = ebfret.analysis.x_lim(x, spread);

    % pick state means to observation mean and variance
    mu = linspace(x_min, x_max, K+2)';
    theta.mu = mu(2:end-1);
    % assume noise proportional to separation of states
    theta.lambda = (12 ./ (x_max - x_min)).^2 * ones(K,1);
    % assume dwell time is 1/2 of average time series length
    theta.tau = 0.5 * mean(T) * ones(K,1);
    % set counts
    counts.mu = 0.025 * strength * ones(K,1);
    counts.lambda = 1 * strength * ones(K,1);
    counts.tau = K * strength * ones(K,1);

    % assign results
    u = ebfret.analysis.hmm.init_prior(theta, counts);
