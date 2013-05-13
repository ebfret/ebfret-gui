function init_analysis(self, num_states)
    if nargin < 2
        num_states = ...
            self.controls.min_states:self.controls.max_states; 
    end
    if ~isprop(self, 'analysis')
        self.analysis = struct('dim', {}, ...
                               'prior', {}, ...
                               'posterior', {}, ...
                               'expect', {}, ...
                               'viterbi', {});
    end
    for k = num_states(num_states>0)
        if (length(self.analysis) < k) || isempty(self.analysis(k).dim)
            self.reset_analysis(k);
        end
    end