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
            % determine locations of states from median
            %
            % we'll calculate succesive median distances away from center
            % and use these values to determine the locations of states,
            % as well as the noise level
            x = self.get_signal();
            x = cat(1, x{:});
            x0 = median(x);
            left(1) = median(x(x<x0));
            right(1) = median(x(x>x0));
            for m = 2:5
                left(m) = median(x(x<left(m-1)));
                right(m) = median(x(x>right(m-1)));
            end

            % get limits from histogram counts
            [h b] = hist(x, linspace(left(end), right(end), min(200, length(x)/10)));
            x_min = max(b(find(h > 0.01 * max(h), 1, 'first')), left(end));
            x_max = min(b(find(h > 0.01 * max(h), 1, 'last')), right(end));

            % pick state means to observation mean and variance
            mu = linspace(x_min, x_max, k+2);
            theta.mu = mu(2:end-1);
            % assume noise level of 0.5 median deviations
            theta.lambda = (4 ./ (right(1) - left(1))).^2 * ones(size(theta.mu));
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