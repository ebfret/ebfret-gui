function clip_outliers(self, x_lim, max_outliers)
    if nargin < 2
        dlg = ...
            inputdlg({'Clip at Minimum:', 'Clip at Maximum:', 'Exclude when number of Outliers exceeds:'}, ...
                'Select Range', 1, {'-0.2', '1.2', '10'});
        x_lim = [str2num(dlg{1}) str2num(dlg{2})];
        max_outliers = str2num(dlg{3});
    end
    for n = 1:length(self.series)
        % clip outlier points
        x = self.series(n).signal;
        self.series(n).signal(x < x_lim(1)) = x_lim(1);
        self.series(n).signal(x > x_lim(2)) = x_lim(2);
        % check if too many outliers
        x = self.series(n).signal(self.series(n).clip.min:self.series(n).clip.max);
        if (sum(x <= x_lim(1)) + sum(x >= x_lim(2))) > max_outliers
            self.series(n).exclude = true;
        else
            (sum(x <= x_lim(1)) + sum(x >= x_lim(2)))
        end
    end
    self.reset_analysis(...
        self.controls.min_states:self.controls.max_states);
    self.set_control(...
        'ensemble', self.controls.ensemble, ...
        'series', self.controls.series);
    self.refresh('ensemble', 'series');
