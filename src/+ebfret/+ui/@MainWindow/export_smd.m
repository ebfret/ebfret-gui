function export_smd(self, varargin)
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
            uiputfile({'*.mat', 'ebFRET SMD export (.mat)'; ...
                       '*.json', 'ebFRET SMD export (.json)'; ...
                       '*.json.gz', 'ebFRET SMD export (.json.gz)';});
        if filebase == 0
            return
        end
        args.filename = sprintf('%s/%s', filepath, filebase);
        switch filetype
            case 1 
                args.filetype = 'mat';
            case 2
                args.filetype = 'json';
            case 3
                args.filetype = 'gz';
        end
    end

    if isempty(args.filetype)
        [void1 void2 args.filetype] = fileparts(args.filename);
    end

    [series, analysis] = self.select_analysis();
    if isempty(series) || isempty(analysis)
        return
    end
    ebfret.io.write_smd(args.filename, series, analysis, 'format', args.filetype);
end