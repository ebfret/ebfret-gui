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
        args.analysis = self.analysis(self.controls.ensemble.min:self.controls.ensemble.max);
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

    for a = 1:length(args.analysis)
        rep = ebfret.analysis.hmm.report(...
                self.get_signal(), ...
                args.analysis(a).prior, ...
                args.analysis(a).expect, ...
                'lowerbound', args.analysis(a).lowerbound, ...
                'splits', args.splits, ...
                'labels', args.labels, ...
                'mapping', args.mapping);
        if a == 1
            report = rep;
        else
            [rep.Series] = deal([]);
            report = cat(1, report(:), rep(:));
        end
    end
    ebfret.io.write_report(args.filename, report);
end