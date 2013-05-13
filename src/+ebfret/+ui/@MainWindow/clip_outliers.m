function clip_outliers(self, x_lim)
    if nargin < 2
        x_lim = ...
            inputdlg({'Clip at Minimum:', 'Clip at Minimum:'}, ...
                'Select Range', 1, {'-0.2', '1.2'});
        x_lim = [str2num(x_lim{1}) str2num(x_lim{2})];
    end
    for n = 1:length(self.series)
        x = self.series(n).signal;
        self.series(n).signal(x < x_lim(1)) = x_lim(1);
        self.series(n).signal(x > x_lim(2)) = x_lim(2);
    end
    self.refresh('ensemble', 'series');
