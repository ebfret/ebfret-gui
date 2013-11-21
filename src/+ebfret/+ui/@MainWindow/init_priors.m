function init_priors(self, theta, counts)
    if nargin < 3
        [theta0, counts0, status] = ebfret.ui.dialog.init_priors
        switch status
            case 0 
                return
            case 1
                K_values = self.controls.min_states:self.controls.max_states;
            case 2
                K_values = self.controls.ensemble.value;
        end
        for k = K_values
            theta(k).mu = linspace(theta0.mu(1), theta0.mu(2), k)';
            theta(k).lambda = ones(k,1) * theta0.lambda;
            theta(k).tau = ones(k,1) * theta0.tau;
            counts(k).mu = ones(k,1) * counts0.mu;
            counts(k).lambda = ones(k,1) * counts0.lambda;
            counts(k).tau = ones(k,1) * counts0.tau;
        end
    end
    for k = 1:length(theta)
        if ~isempty(theta(k).mu)
            self.analysis(k).prior = ebfret.analysis.hmm.init_prior(theta(k), counts(k));
            self.reset_posterior(k);
        end
    end
    self.refresh('ensemble', 'series');
end