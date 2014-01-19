function export_traces(self, varargin)
    ip = inputParser();
    ip.StructExpand = true;
    ip.addOptional('filename', '', @isstr);       
    ip.addParamValue('filetype', '', @isstr);       
    ip.addParamValue('analysis', struct([]), @isstruct);       
    ip.addParamValue('group', '', @isstr);       
    ip.parse(varargin{:});
    args = ip.Results;

    filetype = '';
    if isempty(args.filename)
        [filebase, filepath, filetype] = ...
            uiputfile({'*.dat', 'ebFRET time series (.dat)'; ...
                       '*.mat', 'ebFRET time series (.mat)'});
        if filebase == 0
            return
        end
        args.filename = sprintf('%s/%s', filepath, filebase);
        switch filetype
            case 1 
                args.filetype = 'dat';
            case 2
                args.filetype = 'mat';
        end
    end

    if isempty(args.filetype)
        [void1 void2 args.filetype] = fileparts(args.filename);
    end

    channels = ebfret.ui.dialog.select_channels();
    if ~any(struct2array(channels))
        return
    end

    [series, analysis] = self.select_analysis();
    data = {};
    for n = 1:length(series)
        dat = {};
        if ~series(n).exclude
            range = series(n).crop.min:series(n).crop.max;
            if channels.donor
                dat{end+1} = series(n).donor(range);
            end
            if channels.acceptor
                dat{end+1} = series(n).acceptor(range);
            end
            if channels.fret
               dat{end+1} = series(n).signal(range);
            end
            if channels.viterbi_state
               dat{end+1} = analysis.viterbi(n).state(:);
            end
            if channels.viterbi_mean
               dat{end+1} = analysis.viterbi(n).mean(:);
            end
            data{n} = cat(2, n * ones(length(dat{1}), 1), dat{:});
        end
    end
    traces = cat(1, data{:});
    switch args.filetype
        case 'dat'
            save(args.filename, '-ascii', 'traces');
        case 'mat'
            save(args.filename, 'traces');
    end
end