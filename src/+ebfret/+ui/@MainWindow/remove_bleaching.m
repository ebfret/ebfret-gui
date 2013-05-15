function remove_bleaching(self)
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
    self.reset_analysis(self.controls.min_states:self.controls.max_states);
    self.refresh('ensemble', 'series');
    self.set_control('clip', struct('max', self.series(self.controls.series.value).clip.max));
end