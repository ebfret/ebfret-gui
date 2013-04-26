function refresh_series(self, index)
    status = get(self, 'status');
    plots = get(self, 'plots');
    series = get(self, 'series');
    analysis = get(self, 'analysis');
    handles = get(self, 'handles');
    if length(series) > 0
        if nargin < 2
            n = status.cur_series;
        else 
            n = index;
            status.cur_series = index;
            set_value(handles.seriesControl, index);
            set(self, 'status', status);
        end
        % generate time series plot and observation histogram
        sph = get(handles.seriesPanel, 'handles');
        x = series(n).values;
        t = series(n).index;
        x_lim = get(sph.axes.obs, 'XLim');
        x_bins = linspace(x_lim(1), x_lim(end), ...
                    max(10, min(length(x)/10, 200)));
        if length(analysis) > 0
            g = analysis(status.cur_analysis).expect(n).z;
            z = analysis(status.cur_analysis).viterbi(n).state;
            mu = analysis(status.cur_analysis).posterior(n).mu;
            plots.time = ebfret.plot.time_series(...
                x(:), 'xdata', t(:), ...
                'weights', g, 'state', z(:), 'mean', mu(z(:)), ...
                'color', {plots.colors.obs, ...
                          plots.colors.viterbi, ...
                          plots.colors.state{:}}, ...
                'markersize', 4);
            plots.obs = ebfret.plot.state_obs(...
                x, 'weights', g, ...
                'xdata', x_bins, ...
                'linestyle', '-', ...
                'color', plots.colors.state);
        else
            plots.time = ebfret.plot.time_series(...
                x(:), 'xdata', t(:), ...
                'color', {plots.colors.obs});
            plots.hist = ebfret.plot.state_obs(...
                x, 'xdata', x_bins, ...
                'linestyle', '-', ...
                'color', plots.colors.obs);
        end
        % update time series x-axis limits
        set(sph.axes.time, 'XLim', [min(t), max(t)]);
        % update time series plot
        set_plots(handles.seriesPanel, ...
                  'time', plots.time, ...
                  'obs', plots.obs);
        set(self, 'plots', plots);
        
        % update posterior plots
        fields = {'mean', 'noise', 'dwell'};
        if length(analysis) > 0
            for f = 1:length(fields)
                post = plots.posterior.(fields{f});
                for p = 1:length(post)
                    post(p).ydata = post(p).ydata(:, n);
                end
                set_plots(handles.seriesPanel, ...
                          fields{f}, post);
            end
        else
            clear_plots(handles.seriesPanel, fields);
        end
    else
        clear_plots(handles.seriesPanel, {'time', 'obs'});
    end
end
