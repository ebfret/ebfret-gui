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
            % assign prior
            self.analysis(k).prior = ...
                ebfret.analysis.hmm.guess_prior(self.get_signal(), k);
            % reset posterior
            self.reset_posterior(k);
        end
    end