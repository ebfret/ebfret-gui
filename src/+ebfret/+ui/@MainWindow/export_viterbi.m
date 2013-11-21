function export_viterbi(self, varargin)
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
            uiputfile({'*.dat', 'ebFRET viterbi paths (.dat)'; ...
                       '*.mat', 'ebFRET viterbi paths (.mat)'});
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

    if isempty(args.analysis)
        args.analysis = self.analysis(self.controls.ensemble.value);
    end

    if isempty(args.group)
        ns = 1:length(self.series);
    else
        ns = find(strcmp({self.series.group}, args.group));
    end

    viterbi = [];
    if ~isempty(ns)
        for n = ns
            z = args.analysis.viterbi(n).state; 
            if ~isempty(z)
                m = args.analysis.posterior(n).mu(z);
                viterbi = cat(1, viterbi, cat(2, n*ones(length(z), 1), m, z));
            end
        end
    end

    switch args.filetype
        case 'dat'
            save(args.filename, '-ascii', 'viterbi');
        case 'mat'
            save(args.filename, 'viterbi');
    end
end