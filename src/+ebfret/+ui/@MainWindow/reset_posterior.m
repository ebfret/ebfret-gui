function reset_posterior(self, analysis_index, series_index)
    if nargin < 2
        if self.controls.run_all
            analysis_index = ...
                self.controls.min_states:self.controls.max_states; 
        else
            analysis_index = ...
                self.controls.ensemble.value;
        end
    end
    if nargin < 3
        series_index = ...
            1:length(self.series);
    end
    for a = analysis_index
        for n = series_index
            s = self.series(n);
            if ~s.exclude && (s.clip.max > s.clip.min)
                % reset posterior
                self.analysis(a).posterior(n) = ...
                    self.analysis(a).prior;
                % reset state weights
                self.analysis(a).expect(n).z = ...
                    ones(self.series(n).clip.max-self.series(n).clip.min+1, ...
                         self.analysis(a).dim.states) ./ self.analysis(a).dim.states;
            else
                fields = fieldnames(self.analysis(a).posterior);
                for f = 1:length(fields)
                    self.analysis(a).posterior(n).(fields{f}) = [];
                end
                self.analysis(a).expect(n).z = [];
            end
            % clear viterbi path
            if isfield(self.analysis(a), 'viterbi') 
                self.analysis(a).viterbi(n).state = [];
                self.analysis(a).viterbi(n).mean = [];
            end
       end
    end