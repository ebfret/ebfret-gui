function export_summary(self, varargin)
    ip = inputParser();
    ip.StructExpand = true;
    ip.addOptional('filename', '', @isstr);       
    ip.addOptional('analysis', struct([]), @isstruct);       
    ip.addParamValue('splits', {}, @iscell);       
    ip.addParamValue('labels', {}, @iscell);       
    ip.addParamValue('mapping', [], @isnumeric);       
    ip.parse(varargin{:});
    args = ip.Results;

    if isempty(args.filename)
        [fname, fpath, findex] = ...
            uiputfile({'*.csv', 'ebFRET analysis summary (.csv)';});
        if fname == 0
            return
        end
        args.filename = sprintf('%s/%s', fpath, fname);
    end

    if isempty(args.analysis)
        args.analysis = self.analysis(self.controls.ensemble.value);
    end

    if isempty(args.splits)
        groups = unique({self.series.group});
        if length(groups) > 1
            splits = cellfun(@(g) find(strcmp({self.series.group}, g)), ...
                        groups, 'UniformOutput', false);
        else
            splits = {};
        end
        args.splits = {1:length(self.series), splits{:}};
        if isempty(args.labels)
            args.labels = {'all', groups{:}};
        end
    end

    report = ebfret.analysis.hmm.report(...
                self.get_signal(), ...
                args.analysis.prior, ...
                args.analysis.expect, ...
                'splits', args.splits, ...
                'labels', args.labels, ...
                'mapping', args.mapping);

    ebfret.data.report.to_csv(args.filename, report);
end