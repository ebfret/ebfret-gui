function remove_bleaching(self, method, thresholds)
    if length(self.series) > 0
        if nargin < 2
            [method, thresholds] = ebfret.ui.remove_bleaching_dlg();
        end
        switch method
            case 1
                % remove photobleaching using auto-detected bleaching point
                for n = 1:length(self.series)
                    if ~isempty(self.series(n).donor) && ~isempty(self.series(n).acceptor)
                        id = ebfret.data.fret.photobleach_index(self.series(n).donor);
                        ia = ebfret.data.fret.photobleach_index(self.series(n).acceptor);

                        % sanity check: donor bleaching should result in acceptor bleaching
                        % but we'll allow a few time points tolerance
                        tol = 5;
                        if (ia < (id + tol));
                            self.series(n).clip.max = min(ia, id);
                        else
                            self.series(n).exclude = true;
                        end
                    end
                end
            case 2 
                % remove photobleaching using manual thresholds
                for n = 1:length(self.series)
                    clip_max = length(self.series(n).signal);
                    
                    if ~isnan(thresholds.fret)
                        clip_max = ...
                            min(min(clip_max, ...
                                [find(self.series(n).signal < thresholds.fret, ...
                                    1, 'first'), inf]));
                    end
                    if ~isnan(thresholds.acc)
                        clip_max = ...
                            min(min(clip_max, ...
                                [find(self.series(n).acceptor < thresholds.acc, ...
                                    1, 'first'), inf]));
                    end
                    if ~isnan(thresholds.don)
                        clip_max = ...
                            min(min(clip_max, ...
                                [find(self.series(n).donor < thresholds.don, ...
                                    1, 'first'), inf]));
                    end
                    clip_max
                    if clip_max > self.series(n).clip.min
                        self.series(n).clip.max = clip_max;
                        self.series(n).exclude = false;
                    else
                        self.series(n).exclude = true;
                    end
                end
        end
        self.reset_analysis(self.controls.min_states:self.controls.max_states);
        self.set_control('clip', struct('max', self.series(self.controls.series.value).clip.max));
        self.refresh('ensemble', 'series');
    end
end