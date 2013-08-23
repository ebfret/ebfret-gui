function refresh(self, panel, index)
    % refresh(self, panel, index)
    % 
    % Replots axes in specified plot panel (one of {'series', 'ensemble'})
    import ebfret.analysis.*
    switch panel
        case 'ensemble' 
            controls = get(self, 'controls');
            handles = get(self, 'handles');
            plots = get(self, 'plots');
            series = get(self, 'series');
            analysis = get(self, 'analysis');
            if (length(series) > 0) && (length(analysis) > 0) && any(~[self.series.exclude])
                a = controls.ensemble.value;
                
                % clear plots
                plots = struct();
                
                % set state colors
                self.controls.colors.state = ...
                    ebfret.plot.line_colors(self.analysis(a).dim.states);

                % intialize other colors if unspecified
                if ~isfield(self.controls.colors, 'obs')
                    self.controls.colors.obs = [0.3, 0.3, 0.3];
                end
                if ~isfield(self.controls.colors, 'viterbi')
                    self.controls.colors.viterbi = [0.66, 0.33, 0.33];
                end
                
                % get observations
                x = self.get_signal();

                % get viterbi states
                for n = 1:length(x)
                    if ~self.series(n).exclude && ~isempty(analysis(a).viterbi(n).state)
                        state{n} = analysis(a).viterbi(n).state;
                    else
                        state{n} = zeros(size(x{n}));
                    end
                end

                % concatenate and filter observations and weights
                x = cat(1, x{:});
                state = cat(1, state{:});
                ind = find(isfinite(x));
                x = x(ind);
                state = state(ind,:);

                % get bin positions
                x_bins = ebfret.plot.get_bins(x, 200, min(0.5 ./ sum(~[self.series.exclude]), 1e-2));

                % plot histogram
                plots.ensemble.posterior.obs = ...
                    ebfret.plot.state_obs(...
                        x, 'state', state, ...
                        'xdata', x_bins, ...
                        'num_states', analysis(a).dim.states, ...
                        'color', self.controls.colors.state, ...
                        'linestyle', '-');
                if self.controls.show.prior
                    try
                        % create prior plots
                        u_m = analysis(a).prior.mu;
                        u_beta = analysis(a).prior.beta;
                        u_a = 0.5 .* analysis(a).prior.nu;
                        u_b = 0.5 ./ analysis(a).prior.W;
                        u_alpha = analysis(a).prior.A;
                        try
                            plots.ensemble.prior.mean = ...
                                ebfret.plot.state_mean(u_m, u_beta, u_a, u_b, ...
                                    'xdata', {plots.ensemble.posterior.mean.xdata}, ...
                                    'linestyle', '--', ...
                                    'color', self.controls.colors.state);
                            plots.ensemble.prior.noise = ...
                                ebfret.plot.state_stdev(u_a, u_b, ...
                                    'xdata', {plots.ensemble.posterior.noise.xdata}, ...
                                    'linestyle', '--', ...
                                    'color', self.controls.colors.state);
                            plots.ensemble.prior.dwell = ...
                                ebfret.plot.state_dwell(u_alpha, ...
                                    'xdata', {plots.ensemble.posterior.dwell.xdata}, ...
                                    'linestyle', '--', ...
                                    'color', self.controls.colors.state);
                        catch err
                            plots.ensemble.prior.mean = ...
                                ebfret.plot.state_mean(u_m, u_beta, u_a, u_b, ...
                                    'linestyle', '--', ...
                                    'color', self.controls.colors.state);
                            plots.ensemble.prior.noise = ...
                                ebfret.plot.state_stdev(u_a, u_b, ...
                                    'linestyle', '--', ...
                                    'color', self.controls.colors.state);
                            plots.ensemble.prior.dwell = ...
                                ebfret.plot.state_dwell(u_alpha, ...
                                    'linestyle', '--', ...
                                    'color', self.controls.colors.state);
                        end
                    catch err
                    end
                end

                if self.controls.show.posterior
                    % create posterior plots
                    try
                        % update all posterior plots
                        ns = find(~[self.series.exclude]);
                        w_m = cat(2, analysis(a).posterior(ns).mu);
                        w_beta = cat(2, analysis(a).posterior(ns).beta);
                        w_a = 0.5 * cat(2, analysis(a).posterior(ns).nu);
                        w_b = 0.5 ./ cat(2, analysis(a).posterior(ns).W);
                        w_alpha = cat(3, analysis(a).posterior(ns).A);
                        plots.ensemble.posterior.mean = ...
                            ebfret.plot.mean(...
                                ebfret.plot.state_mean(w_m, w_beta, w_a, w_b, ...
                                                       'linestyle', '-', ...
                                                       'color', self.controls.colors.state));
                        plots.ensemble.posterior.noise = ...
                            ebfret.plot.mean(...
                                ebfret.plot.state_stdev(w_a, w_b, ...
                                                        'linestyle', '-', ...
                                                        'color', self.controls.colors.state));
                        E_tau = -1 ./ log(diag(ebfret.normalize(mean(w_alpha,3), 2)));
                        plots.ensemble.posterior.dwell = ...
                            ebfret.plot.mean(...
                                ebfret.plot.state_dwell(w_alpha, ...
                                                        'xdata', arrayfun(@(tau) exp(linspace(log(0.01 * tau), log(100 * tau), 101))', ...
                                                                    E_tau(:)', 'UniformOutput', false), ...
                                                        'linestyle', '-', ...
                                                        'color', self.controls.colors.state));
                    catch err
                        rethrow(err)
                    end
                end

                % store prior and posterior plots
                set(self, 'plots', plots);

                % set axis properties
                sph = get(handles.seriesPanel, 'handles');
                eph = get(handles.ensemblePanel, 'handles');

                % set logarithmic axis for dwell time plots
                set(eph.axes.dwell, 'XScale', 'log');

                % set axis limits for observations
                [x_lim, y_lim] = ...
                    ebfret.plot.get_lim(plots.ensemble.posterior.obs, ...
                        1e-2, [0.05, 0.05, 0.05, 0.15]);
                try
                    set(sph.axes.signal, 'YLim', x_lim);
                    % set(sph.axes.obs, ...
                    %     'XLim', x_lim, 'YLim', y_lim, ...
                    %     'YTick', linspace(y_lim(1), y_lim(2), 5));
                    set(eph.axes.obs, ...
                        'XLim', x_lim, 'YLim', y_lim);
                catch
                    % copy limits from histogram plot
                    x_lim = get(eph.axes.obs, 'xlim');
                    set(sph.axes.signal, 'YLim', x_lim);
                end
                   
                % update ensemble plots
                handles = get(self, 'handles');
                set_plots(handles.ensemblePanel, ...
                    'obs', plots.ensemble.posterior.obs);
                
                % get occupancy dependent scaling
                if self.controls.scale_plots
                    scale = mean(cat(2, analysis(a).expect.z),2);
                    if isempty(scale)
                        scale = 1;
                    end
                else
                    scale = 1;
                end

                axes = {'mean', 'noise', 'dwell'};
                thresholds = [0.02, 0.02, 0.001];
                for ax = 1:length(axes)
                    clear_plots(handles.ensemblePanel, axes{ax});
                    ax_plots = struct([]);
                    try 
                        ax_plots = cat(1, ax_plots, ...
                            ebfret.plot.scale(plots.ensemble.prior.(axes{ax})(:), scale));
                    catch
                    end
                    try 
                        ax_plots = cat(1, ax_plots, ...
                            ebfret.plot.scale(plots.ensemble.posterior.(axes{ax})(:), scale));
                    catch
                    end
                    if ~isempty(ax_plots)
                        % set axis limits
                        [x_lim, y_lim] = ...
                            ebfret.plot.get_lim(ax_plots, ...
                                thresholds(ax), [0.05, 0.05, 0.05, 0.25]);
                        if strcmp(axes{ax}, 'dwell')
                            x_lim = [max(1e-4 * x_lim(2), x_lim(1)), x_lim(2)];
                            % y_lim = [2e-3 * y_lim(2), 2 * y_lim(2)];
                        end
                        if all(isfinite(x_lim))
                            set(eph.axes.(axes{ax}), 'XLim', x_lim);
                        end
                        if all(isfinite(y_lim))
                            set(eph.axes.(axes{ax}), 'YLim', y_lim);
                        end
                        % update plots
                        set_plots(handles.ensemblePanel, ...
                            axes{ax}, ax_plots);
                    end
                end

                axes = struct2cell(eph.axes);
                axes = [axes{:}];
                for ax = axes
                    % x-tick label formatting
                    xtick = get(ax, 'xtick');
                    set(ax, 'xticklabel', ebfret.num_to_str(xtick));
                    % y-tick label formatting
                    ylim = get(ax, 'ylim');
                    set(ax, 'YTick', linspace(ylim(1), ylim(2), 5));
                    set(ax, 'YTickLabel', {});
                end

                % update time series plots
                self.refresh('series');
            else
                clear_plots(handles.ensemblePanel, ...
                    {'obs', 'mean', 'noise', 'dwell'});
            end
        case 'series'
            controls = get(self, 'controls');
            plots = get(self, 'plots');
            series = get(self, 'series');
            handles = get(self, 'handles');
            if length(self.analysis) > 0
                analysis = self.analysis(self.controls.ensemble.value);
            else
                analysis = struct([]);
            end
            
            clear_plots(handles.seriesPanel);
            if (length(series) > 0) & (self.controls.series.value > 0)
                n = self.controls.series.value;
                plots.series = struct();
                sph = get(handles.seriesPanel, 'handles');
                % get crop range
                if ~self.series(n).exclude
                    crop = self.series(n).crop;
                else
                    crop.max = length(self.series(n).signal);
                    crop.min = crop.max;
                end
                pargs = struct('crop', self.series(n).crop, 'markersize', 4);
                pargs.color = {self.controls.colors.obs}; 
                if self.controls.show.viterbi
                    try
                        % set state colors
                        self.controls.colors.state = ...
                            ebfret.plot.line_colors(analysis.dim.states);
                        % use viterbi path if available
                        pargs.state = analysis.viterbi(n).state;
                        pargs.num_states = analysis.dim.states;
                        pargs.color = {self.controls.colors.obs, ...
                                       self.controls.colors.viterbi, ...
                                       self.controls.colors.state{:}};
                    catch
                        % just plot signal otherwise
                    end
                end
                % plot signal
                plots.series.signal = ebfret.plot.time_series(...
                    series(n).signal, series(n).time, pargs);
                try
                    % plot raw signal if available
                    pargs.color{1} = self.controls.colors.donor;
                    donor = ebfret.plot.time_series(...
                                series(n).donor, series(n).time, pargs);
                    pargs.color{1} = self.controls.colors.acceptor;
                    acceptor = ebfret.plot.time_series(...
                                    series(n).acceptor, series(n).time, pargs);
                    if self.controls.show.viterbi
                        plots.series.raw = [donor(1), donor(3:end), acceptor(1), acceptor(3:end)];
                    else
                        plots.series.raw = [donor(:)', acceptor(:)'];
                    end
                catch err
                    % simply ignore errors
                    rethrow(err)
                end
                
                % update time series x-axis limits
                t_min = min(self.series(n).time);
                t_max = self.series(n).time(min(length(self.series(n).signal), ...
                                                    crop.max + self.controls.crop_margin));
                set(sph.axes.signal, 'XLim', [t_min, t_max]);
                set(sph.axes.raw, 'XLim', [t_min, t_max]);

                % tick label formatting
                axes = struct2cell(sph.axes);
                axes = [axes{:}];
                for ax = axes
                    xtick = get(ax, 'xtick');
                    set(ax, 'xticklabel', ebfret.num_to_str(xtick));
                end

                % update time series plot
                set_plots(handles.seriesPanel, plots.series)
                set(self, 'plots', plots);
            else
                clear_plots(handles.seriesPanel, ...
                    {'signal', 'raw'});
            end
    end
end