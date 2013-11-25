function clip_outliers(self, x_lim, max_outliers)
    if nargin < 2
        dlg = ...
            inputdlg({'Clip at Minimum:', 'Clip at Maximum:', 'Exclude when number of Outliers exceeds:'}, ...
                'Select Range', 1, {'-0.2', '1.2', '10'});
        % return if user hit cancel
        if isempty(dlg)
            return
        end
        x_lim = [str2num(dlg{1}) str2num(dlg{2})];
        max_outliers = str2num(dlg{3});
    end
    total_points = 0;
    total_out = 0;
    removed = 0;
    for n = 1:length(self.series)
        % set clipping controls
        self.set_control('clip', struct('min', x_lim(1), 'max', x_lim(2)));
        % check if too many outliers
        x = self.series(n).signal(self.series(n).crop.min:self.series(n).crop.max);
        out = sum(x <= x_lim(1)) + sum(x >= x_lim(2));

        if out > max_outliers
            self.series(n).exclude = true;
            removed = removed + 1;
        end
        total_points = total_points + length(x);
        total_out = total_out + out;
    end
    if nargin < 2
        msgbox({sprintf('Found %d outlier points (%.3f%% of total).', total_out, 100 * total_out / total_points); ...
                ''; ...
                sprintf('Excluded %d out of %d time series from analysis.', removed, length(self.series))});
    end
    % reset posteriors
    self.reset_posterior(self.controls.min_states:self.controls.max_states);
    % ask to update priors
    self.update_priors()
    % update plots
    self.refresh('ensemble', 'series');
