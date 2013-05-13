function save_data(self)
    [fname, fpath, findex] = ...
        uiputfile({'*.mat', 'ebFRET saved session (.mat)';});
    switch findex
        case 1
            controls = get(self, 'controls');
            series = get(self, 'series');
            analysis = get(self, 'analysis');
            plots = get(self, 'plots');
            save(sprintf('%s/%s', fpath, fname), ...
                'controls', 'series', 'analysis', 'plots');
    end
end
