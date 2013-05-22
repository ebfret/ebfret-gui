function load_data(self, file, ftype)
    if nargin < 2
        [fname, fpath, ftype] = ...
            uigetfile({'*.mat', 'ebFRET saved session (.mat)';
                       '*.dat', 'Raw donor-acceptor time series (.dat)';
                       '*.tsv', 'SF-Tracer donor-acceptor time series (.tsv)';});
        file = sprintf('%s/%s', fpath, fname);
    end
        % uigetfile({'*.mat', 'ebFRET saved session (.mat)';, ...
        %            '*.mat', 'vbFRET saved session (.mat)';, ... 
        %            '*.mat;*.smd', 'Single-molecule Data Format (.smd,.mat)';, ... 
        %            '*.tsv', 'SF-Tracer donor-acceptor time series (.tsv)'; ...
        %            '*.dat', 'Raw donor-acceptor time series (.dat)'});
    if ~isempty(self.series)
        choice = questdlg('Would you like to keep already loaded time series, or replace them?\n\n Warinng: Previous analysis will be lost.', ...
                    'Append Data?', 'Keep', 'Replace', 'Cancel', 'Keep');
        switch choice
            case 'Keep'
                append = true;
            case 'Replace'
                append = false;
            case 'Cancel'
                return
        end
    else
        append = false;
    end
    switch ftype
        case 1
            session = load(file);
            self.series = session.series; 
            self.analysis = session.analysis;
            self.plots = session.plots;
            self.set_control(...
                'ensemble', session.controls.ensemble, ...
                'series', session.controls.series);
            self.set_control(rmfield(session.controls, {'ensemble', 'series'}));
        case 2
            [don acc labels] = ...
                ebfret.data.fret.load_raw(file, 'has_labels', true);
        case 3
            [don acc] = ...
                ebfret.data.fret.load_sf_tracer(file);
            labels = 1:length(don);
    end

    if (ftype == 2) || (ftype == 3)
        % strip empty time series
        ns = find(~(cellfun(@isempty, don)) & ~(cellfun(@isempty, acc)));
        don = {don{ns}};
        acc = {acc{ns}};
        labels = labels(ns);

        % format labels as string
        labels = arrayfun(...
                    @(l) sprintf('%d', l), labels, ...
                    'UniformOutput', false);

        % initialize series struct
        series = struct('file', {}, ...
                        'label', {}, ...
                        'time', {}, ...
                        'signal', {}, ...
                        'donor', {}, ...
                        'acceptor', {}, ...
                        'clip', {}, ...
                        'exclude', {});
        
        for n = 1:length(don)
            [void series(n).file] = fileparts(file);
            series(n).label = labels{n};
            series(n).donor = don{n}(:);
            series(n).acceptor = acc{n}(:);
            series(n).time = (1:length(series(n).donor))';
            series(n).signal = series(n).acceptor ...
                                ./ (series(n).acceptor + series(n).donor);
            series(n).exclude = false;
            series(n).clip.min = 1;
            series(n).clip.max = length(series(n).time);
        end
        if ~append
            self.series = series;
        else
            self.series =cat(1, self.series(:), series(:));
        end
        self.reset_analysis(self.controls.min_states:self.controls.max_states);
        self.set_control(...
            'series', struct(...
                        'min', 1, ...
                        'max', length(self.series), ...
                        'value', 1));
        self.set_control('ensemble', struct('value', self.controls.min_states));
        self.set_control('run_analysis', false);
    end
    self.refresh('ensemble', 'series');
end
