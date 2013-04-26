function save_data(self)
    [fname, fpath, findex] = ...
        uiputfile({'*.mat', 'ebFRET saved session (.mat)';});
    switch findex
        case 1
            status = get(self, 'status');
            series = get(self, 'series');
            analysis = get(self, 'analysis');
            plots = get(self, 'plots');
            save(sprintf('%s/%s', fpath, fname), ...
                'status', 'series', 'analysis', 'plots');
    end
end
