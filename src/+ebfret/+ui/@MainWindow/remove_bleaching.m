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

                        % % sanity check: donor bleaching should result in acceptor bleaching
                        % % but we'll allow a few time points tolerance
                        % tol = 5;
                        % if (ia < (id + tol));
                        %     self.series(n).clip.max = min(ia, id);
                        % else
                           %     self.series(n).exclude = true;
                        % end
                    end
                end
            case 2 
                % remove photobleaching using manual thresholds
                W = 7;
                W0 = floor(W/2);
                kernel = ebfret.analysis.normalize(exp(-linspace(-1.5,1.5,W).^2));
                for n = 1:length(self.series)
                    clip_max = length(self.series(n).signal) - self.series(n).clip.min;
                    if ~isnan(thresholds.fret)
                        signal = conv(self.series(n).signal(self.series(n).clip.min:end), kernel, 'valid');
                        clip_max = ...
                            min(min(clip_max, ...
                                [find(signal < thresholds.fret, 1, 'first') + W0, inf]));
                    end
                    if ~isnan(thresholds.acc) || ~isnan(thresholds.sum)
                        acceptor = conv(self.series(n).acceptor(self.series(n).clip.min:end), kernel, 'valid');
                        clip_max = ...
                            min(min(clip_max, ...
                                [find(acceptor < thresholds.acc, 1, 'first') + W0, inf]));
                    end
                    if ~isnan(thresholds.don) || ~isnan(thresholds.sum)
                        donor = conv(self.series(n).donor(self.series(n).clip.min:end), kernel, 'valid');
                        clip_max = ...
                            min(min(clip_max, ...
                                [find(donor < thresholds.don, 1, 'first') + W0, inf]));
                    end
                    if ~isnan(thresholds.sum)
                        clip_max = ...
                            min(min(clip_max, ...
                                [find(donor+acceptor < thresholds.sum, 1, 'first') + W0, inf]));
                    end
                    if ~isnan(thresholds.pad)
                        clip_max = ...
                            clip_max - thresholds.pad;
                    end
                    if clip_max > 0
                        self.series(n).clip.max = clip_max + self.series(n).clip.min;
                        self.series(n).exclude = false;
                    else
                        self.series(n).exclude = true;
                    end
                end
            otherwise
                return
        end
        self.reset_analysis(self.controls.min_states:self.controls.max_states);
        self.set_control('clip', struct('max', self.series(self.controls.series.value).clip.max));
        self.refresh('ensemble', 'series');
    end
end