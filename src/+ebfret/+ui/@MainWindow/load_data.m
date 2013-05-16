function load_data(self)
    [fname, fpath, ftype] = ...
        uigetfile({'*.mat', 'ebFRET saved session (.mat)';
                   '*.dat', 'Raw donor-acceptor time series (.dat)';
                   '*.tsv', 'SF-Tracer donor-acceptor time series (.tsv)';});
        % uigetfile({'*.mat', 'ebFRET saved session (.mat)';, ...
        %            '*.mat', 'vbFRET saved session (.mat)';, ... 
        %            '*.mat;*.smd', 'Single-molecule Data Format (.smd,.mat)';, ... 
        %            '*.tsv', 'SF-Tracer donor-acceptor time series (.tsv)'; ...
        %            '*.dat', 'Raw donor-acceptor time series (.dat)'});
    switch ftype
        case 1
            session = load(sprintf('%s/%s', fpath, fname));
            self.series = session.series; 
            self.analysis = session.analysis;
            self.plots = session.plots;
            self.set_control(...
                'ensemble', session.controls.ensemble, ...
                'series', session.controls.series);
            self.set_control(rmfield(session.controls, {'ensemble', 'series'}));
        case 2
            [don acc] = ebfret.data.fret.load_raw(...
                            sprintf('%s/%s', fpath, fname), ...
                            'has_labels', true);
        case 3
            [don acc] = ebfret.data.fret.load_sf_tracer(...
                            sprintf('%s/%s', fpath, fname));
    end

    if (ftype == 2) || (ftype == 3)


        % [don acc] = ebfret.data.fret.remove_bleaching(...
        %                 'donor', don, 'acceptor', acc);
        % don = don(~cellfun(@isempty, don));
        % acc = acc(~cellfun(@isempty, acc));
        % for n = 1:length(data)
        %     data(n).fret = (data(n).acceptor+eps) ./ (data(n).donor + data(n).acceptor + eps);
        %     data(n).fret(data(n).fret>1.5) = 1.5; 
        %     data(n).fret(data(n).fret<-0.5) = -0.5; 
        % end 
        self.series = struct('time', {}, ...
                             'signal', {}, ...
                             'donor', {}, ...
                             'acceptor', {}, ...
                             'clip', {}, ...
                             'exclude', {});
        for n = 1:length(don)
            self.series(n).donor = don{n}(:);
            self.series(n).acceptor = acc{n}(:);
            self.series(n).time = (1:length(self.series(n).donor))';
            self.series(n).signal = self.series(n).acceptor ...
                                ./ (self.series(n).acceptor + self.series(n).donor);
            % self.series(n).signal(self.series(n).signal<-0.2) = -0.2;                                
            % self.series(n).signal(self.series(n).signal>1.2) = 1.2;                                
            self.series(n).exclude = false;
            self.series(n).clip.min = 1;
            self.series(n).clip.max = length(self.series(n).time);
        end

        self.reset_analysis(self.controls.min_states:self.controls.max_states);
        self.set_control('ensemble', struct('value', self.controls.min_states));
        self.set_control(...
            'series', struct(...
                        'min', 1, ...
                        'max', length(self.series), ...
                        'value', 1));
        self.set_control('run_analysis', false);
    end
    self.refresh('ensemble', 'series');
end
