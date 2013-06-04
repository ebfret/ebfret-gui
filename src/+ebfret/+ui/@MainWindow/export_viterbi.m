function export_viterbi(self, varargin)
    ip = inputParser();
    ip.StructExpand = true;
    ip.addOptional('filename', '', @isstr);       
    ip.addOptional('analysis', struct([]), @isstruct);       
    ip.addOptional('group', '', @isstr);       
    ip.parse(varargin{:});
    args = ip.Results;

    if isempty(args.filename)
        [fname, fpath, findex] = ...
            uiputfile({'*.dat', 'ebFRET viterbi paths (.dat)';});
        if fname == 0
            return
        end
        args.filename = sprintf('%s/%s', fpath, fname);
    end

    if isempty(args.analysis)
        args.analysis = self.analysis(self.controls.ensemble.value);
    end

    if isempty(args.group)
        ns = 1:length(self.series);
    else
        ns = find(strcmp({self.series.group}, args.group));
    end

    data = [];
    if ~isempty(ns)
        for n = ns
            z = args.analysis.viterbi(n).state; 
            if ~isempty(z)
                m = args.analysis.posterior(n).mu(z);
                data = cat(1, data, cat(2, n*ones(length(z), 1), m, z));
            end
        end
    end

    save(args.filename, '-ascii', 'data');
end