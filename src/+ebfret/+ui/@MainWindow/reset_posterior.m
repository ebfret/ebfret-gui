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
            % reset posterior
            self.analysis(a).posterior(n) = ...
                self.analysis(a).prior;
            self.analysis(a).expect(n) = ...
                struct('z', [], 'z1', [], 'zz', [], 'x', [], 'xx', []);
            if s.exclude || (s.crop.max == s.crop.min)
                fields = fieldnames(self.analysis(a).posterior);
                for f = 1:length(fields)
                    self.analysis(a).posterior(n).(fields{f}) = [];
                end
            end
            % clear lower bound
            self.analysis(a).lowerbound(n) = 0;
            % clear viterbi path
            self.analysis(a).viterbi(n).state = [];
            self.analysis(a).viterbi(n).mean = [];
       end
    end