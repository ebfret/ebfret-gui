function x = get_signal(self, series_index)
    if nargin < 2
        series_index = 1:length(self.series);
    end
    x = {};
    for n = series_index
        s = self.series(n);
        if ~s.exclude
            x{n} = s.signal(s.clip.min:s.clip.max);
            % % remove inf and nan instances
            % x{n}(x{n} == inf) = 1 ./ eps;
            % x{n}(x{n} == -inf) = -1 ./ eps;
            % x{n}(isnan(x{n})) = 0;
        end
    end
    if isscalar(series_index) 
        if ~isempty(x)
            x = x{n};
        else
            x = [];
        end
    end