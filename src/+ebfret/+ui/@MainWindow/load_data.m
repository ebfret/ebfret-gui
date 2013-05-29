function load_data(self, files, ftype)
    if nargin < 2
        [fname, fpath, ftype] = ...
            uigetfile({'*.mat', 'ebFRET saved session (.mat)';
                       '*.dat', 'Raw donor-acceptor time series (.dat)';
                       '*.tsv', 'SF-Tracer donor-acceptor time series (.tsv)';}, ...
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
        self.series = session.series; 
        self.analysis = session.analysis;
        self.plots = session.plots;
        self.set_control(...
            'ensemble', session.controls.ensemble, ...
            'series', session.controls.series);
        self.set_control(rmfield(session.controls, ...
            {'ensemble', 'series', 'run_analysis'}));
    end

    if (ftype == 2) || (ftype == 3)
        dlg = waitbar(0,'Loading','Name','Loading datasets', ...
                      'CreateCancelBtn', 'cancelled = true', ...
                      'Interpreter', 'none');
        % set(dlg,'DefaultTextInterpreter','none');
        for f = 1:length(files)
            fname = fileparts(files{f})
            waitbar((f-1)/(length(files)-1), dlg, sprintf('Loading: %s', fname));
            if (ftype == 2)
                [dons{f} accs{f} labels{f}] = ...
                    ebfret.data.fret.load_raw(files{f}, 'has_labels', true);
            elseif (ftype == 3)
                [dons{f} accs{f}] = ...
                    ebfret.data.fret.load_sf_tracer(files{f});
                labels{f} = 1:length(dons{f});
            end
            if cancelled
                return
            end
        end

        waitbar((f-1)/(length(files)-1), dlg, ...
            ebfret.escape_tex(...
                sprintf('Read %d time series from %d files', ...
                    sum(cellfun(@length, labels)), length(files))));

        if ~append
            self.series = struct([]);
            group = 1;
        else
            group = max([self.series.group]) + 1;
        end
        for f = 1:length(files)
            % strip empty time series
            ns = find(~(cellfun(@isempty, dons{f})) & ~(cellfun(@isempty, accs{f})));
            don = {dons{f}{ns}};
            acc = {accs{f}{ns}};
            label = labels{f}(ns);

            % format labels as string
            label = arrayfun(...
                        @(l) sprintf('%d', l), label, ...
                        'UniformOutput', false);

            % initialize series struct
            series = struct('file', {}, ...
                            'label', {}, ...
                            'group', {}, ...
                            'time', {}, ...
                            'signal', {}, ...
                            'donor', {}, ...
                            'acceptor', {}, ...
                            'clip', {}, ...
                            'exclude', {});
            
            for n = 1:length(don)
                [void series(n).file] = fileparts(files{f});
                series(n).label = label{n};
                series(n).group = group;
                series(n).donor = don{n}(:);
                series(n).acceptor = acc{n}(:);
                series(n).time = (1:length(series(n).donor))';
                series(n).signal = series(n).acceptor ...
                                    ./ (series(n).acceptor + series(n).donor);
                series(n).exclude = false;
                series(n).clip.min = 1;
                series(n).clip.max = length(series(n).time);
            end
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
        delete(dlg);
    end
    self.set_control('run_analysis', false);
    self.refresh('ensemble', 'series');
end
