function u = guess_prior(x, K, quantiles, strength)
	% u = guess_prior(x, K, quantiles)
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
    % quantiles : (1 x 2) (default: [0.01, 0.99])
    %   Quantiles for positioning min and max state means
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
        quantiles = [0.01, 0.99];
    end
    if nargin < 4
        strength = 1;
    end
    T = cellfun(@length, x);
    x = cat(1, x{:});
    x_min = quantile(x, quantiles(1));
    x_max = quantile(x, quantiles(2));

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
