function update_priors(self, choice)
    if nargin < 2
        choice =  questdlg('Update prior distributions ?', ...
                           'Update Priors', ...
                           'Auto', 'Manual', 'Keep Current', 'Keep Current');
    end
    switch lower(choice)
        case 'auto'
            num_states = ...
                self.controls.min_states:self.controls.max_states; 
            for k = num_states(num_states>0)
                self.analysis(k).prior = ...
                    ebfret.analysis.hmm.guess_prior(self.get_signal(), k);
            end
        case 'manual'
            self.init_priors()
    end
