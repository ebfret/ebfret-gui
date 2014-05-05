function load_data(self, files, ftype)
    if nargin < 2
        [fname, fpath, ftype] = ...
            uigetfile({'*.mat', 'ebFRET saved session (.mat)';
                       '*.dat', 'Raw donor-acceptor time series (.dat)';
                       '*.tsv', 'SF-Tracer donor-acceptor time series (.tsv)';
                       '*.mat', 'SMD time series (.mat)';
                       '*.json', 'SMD time series (.json)';
                       '*.json.gz', 'SMD time series (.json.gz)';}, ...
                       'multiselect', 'on');
        if (ftype == 0)
            return
        end
        if iscell(fname)
            for f = 1:length(fname)
                files{f} = sprintf('%s/%s', fpath, fname{f});
            end
        else
            files{1} = sprintf('%s/%s', fpath, fname);
        end
    end
    if isstr(files)
        files = {files};
    end
    % uigetfile({'*.mat', 'ebFRET saved session (.mat)';, ...
    %            '*.mat', 'vbFRET saved session (.mat)';, ... 
    %            '*.mat;*.smd', 'Single-molecule Data Format (.smd,.mat)';, ... 
    %            '*.tsv', 'SF-Tracer donor-acceptor time series (.tsv)'; ...
    %            '*.dat', 'Raw donor-acceptor time series (.dat)'});
    if ~isempty(self.series) & (ftype ~= 1)
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
    % this seems a bit of a cludge, but works
    cancelled = false;
    function status = cancel()
        cancelled = true; 
        status = true;
    end
    if ftype == 1
        if length(files) > 1
            warndlg('You cannot load more than one ebfret session at a time. The first selected file will be loaded.');
        end
        session = load(files{1});
        % check for old style group label (numeric instead of string)
        if all(cellfun(@isnumeric, {session.series.group}))
            groups = unique([session.series.group]);
            ns = arrayfun(@(g) find([session.series.group] == g), ...
                    groups, 'UniformOutput', false);
            for g = 1:length(groups)
                [session.series(ns{g}).group] = deal(sprintf('group %d', groups(g)));
            end
        end
        self.series = session.series; 
        self.analysis = session.analysis;
        % self.plots = session.plots;
        self.set_control(...
            'ensemble', session.controls.ensemble, ...
            'series', session.controls.series);
        self.set_control(rmfield(session.controls, ...
            {'ensemble', 'series', 'run_analysis'}));
    end

    if any(ftype == 2:6)
        dlg = waitbar(0, 'Loading', 'Name', 'Loading datasets', ...
                         'CreateCancelBtn', 'cancelled = true', ...
                         'Interpreter', 'none');
        try 
            for f = 1:length(files)
                [void name] = fileparts(files{f});
                waitbar((f-1)/(length(files)-1), dlg, ebfret.escape_tex(name));
                switch ftype
                case 2
                    try
                        [dons{f} accs{f} labels{f}] = ...
                            ebfret.io.load_raw(files{f}, 'has_labels', true);
                    catch err
                        error('ebfret:load_data:wrong_format', ...
                              'File "%s" could not be loaded as Raw data. Type "help ebfret.data.load_raw" for a description of supported formats.', files{f});
                    end
                case 3
                    try
                        [dons{f} accs{f}] = ...
                            ebfret.io.load_sf_tracer(files{f});
                        labels{f} = 1:length(dons{f});
                    catch err
                        error('ebfret:load_data:wrong_format', ...
                              'File "%s" could not be loaded as SF-Tracer data.', files{f});
                    end
                case 4
                    smd{f} = load(files{f});
                case 5
                    smd{f} = ebfret.io.load_json(files{f});
                case 6
                    smd{f} = ebfret.io.load_json(files{f}, 'gzip', true);
                end
                if cancelled
                    return
                end
            end

            switch ftype
            case {2,3}
                num_series = sum(cellfun(@length, labels));
            case {4,5,6}
                num_series = sum(cellfun(@(s) length(s.data), smd));
            end
            waitbar((f-1)/(length(files)-1), dlg, ...
                ebfret.escape_tex(...
                    sprintf('Read %d time series from %d files', ...
                        num_series, length(files))));

            if ~append
                self.series = struct([]);
                group = 'group 1';
            else
                num_groups = length(unique({self.series.group}));
                group = sprintf('group %d', num_groups+1);
            end
            for f = 1:length(files)
                % initialize series struct
                series = struct('file', {}, ...
                                'label', {}, ...
                                'group', {}, ...
                                'time', {}, ...
                                'signal', {}, ...
                                'donor', {}, ...
                                'acceptor', {}, ...
                                'crop', {}, ...
                                'exclude', {});

                switch ftype
                case {2,3}
                    % strip empty time series
                    ns = find(~(cellfun(@isempty, dons{f})) & ~(cellfun(@isempty, accs{f})));
                    don = {dons{f}{ns}};
                    acc = {accs{f}{ns}};
                    label = labels{f}(ns);

                    % format labels as string
                    label = arrayfun(...
                                @(l) sprintf('%d', l), label, ...
                                'UniformOutput', false);

                    
                    for n = 1:length(don)
                        [void series(n).file] = fileparts(files{f});
                        series(n).label = label{n};
                        series(n).group = group;
                        series(n).donor = don{n}(:);
                        series(n).acceptor = acc{n}(:);
                        series(n).time = (1:length(series(n).donor))';
                        series(n).signal = (series(n).acceptor + eps) ...
                                            ./ (series(n).acceptor + series(n).donor + eps);
                        series(n).crop.min = 1;
                        series(n).crop.max = length(series(n).time);
                    end
                case {4,5,6}
                    channels = ebfret.ui.dialog.assign_smd_channels(smd{f}.columns);
                    series = struct([]);
                    data = smd{f}.data
                    for n = 1:length(data)
                        [void series(n).file] = fileparts(files{f});
                        series(n).label = data(n).id;
                        series(n).time = data(n).index;
                        if channels.fret == channels.fret
                            series(n).donor = zeros(size(data(n).index));
                            series(n).acceptor = zeros(size(data(n).index));
                            series(n).signal = data(n).values(:, channels.fret);
                        else
                            series(n).donor = data(n).values(:, channels.donor);
                            series(n).acceptor = data(n).values(:, channels.acceptor);
                            series(n).signal = (series(n).acceptor + eps) ./ (series(n).acceptor + series(n).donor + eps);
                        end
                        series(n).crop.min = 1;
                        series(n).crop.max = length(series(n).time);
                    end
                end
                [series.exclude] = deal(false);
                if isempty(self.series)
                    self.series = series;
                else
                    self.series = cat(1, self.series(:), series(:));
                end
            end

            self.reset_analysis(self.controls.min_states:self.controls.max_states);
            self.set_control(...
                'series', struct(...
                            'min', 1, ...
                            'max', length(self.series), ...
                            'value', 1));
            self.set_control('ensemble', struct('value', self.controls.min_states));
        catch err
            err_dlg = errordlg(err.message);
        end
        delete(dlg);
    end

    % this is to ensure analysis does not start immediately
    self.set_control('run_analysis', false);
    % this updates enabled/disabled status of the controls in the view menu
    self.set_control('show', struct());
    self.set_control('scale_plots', self.controls.scale_plots);
    % this does a replot for good measure
    self.refresh('ensemble', 'series');
end
