function reset_analysis(self, num_states)
    if nargin < 2
        if self.controls.run_all
            num_states = ...
                self.controls.min_states:self.controls.max_states; 
        else
            num_states = ...
                self.controls.ensemble.value;
        end
    end
    for k = num_states(num_states>0)
        % set number of states
        self.analysis(k).dim.states = k;
        % initialize prior and posterior (if we have time series data)
        if (length(self.series) > 0) 
            % get mean and variance of observations
            x = self.get_signal();
            x = cat(1, x{:});
            idx = find(isfinite(x));
            x = x(idx);
            % mean and variance of observations
            E_x = mean(x);
            V_x = var(x);
            % pick state means to observation mean and variance
            mu = sqrt(V_x) * (1:k)';
            theta.mu = mu - mean(mu) + E_x;
            % assume noise level of 0.5 standard deviations
            theta.lambda = 4 ./ V_x .* ones(size(theta.mu));
            % assume dwell time is 1/10th of average time series length
            theta.tau = 0.1 * mean(cellfun(@length, {self.series.signal})) ...
                                   * ones(size(theta.mu));
            % set counts
            counts.mu = 0.1 * ones(size(theta.mu));
            counts.lambda = 1 * ones(size(theta.lambda));
            counts.tau = length(theta.tau) * ones(size(theta.tau));
            % assign prior
            self.analysis(k).prior = ...
                ebfret.analysis.hmm.init_prior(theta, counts);
            % reset posterior
            self.reset_posterior(k);
        end
    end