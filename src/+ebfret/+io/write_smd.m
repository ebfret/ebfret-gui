function write_smd(filename, series, analysis, varargin)
    % write_smd(filename, series, analysis, varargin)
    % 
    % Writes time series and analysis to single-molecule dataset (SMD).
    %
    % Inputs
    % ------
    %
    %   filename : string
    %       Output file
    %
    %   series : struct
    %       Time series data (as stored in ebfret.ui.MainWindow)
    %
    %   analysis : struct
    %       Analysis data (as stored in ebfret.ui.MainWindow)
    %
    %
    % Variable Inputs
    % ---------------
    %
    %   format : {'mat', 'json', 'json.gz'}
    %       File format to use. Guessed from file name if not specified. 
    ip = inputParser();
    ip.StructExpand = true;
    ip.KeepUnmatched = true;
    ip.addRequired('filename', @isstr);
    ip.addRequired('series', @isstruct);
    ip.addRequired('analysis', @isstruct);
    ip.addParamValue('format', '',  @(s) any(strcmpi(s, {'', 'mat', 'json', 'gz'})));
    ip.parse(filename, series, analysis, varargin{:});
    args = ip.Results;
    opts = ip.Unmatched;

    % determine output file format
    if isempty(args.format)
        [p, f, args.format] = fileparts(args.filename);
        args.format = args.format(2:end);
        if ~any(strcmpi(args.format, {'mat', 'json', 'gz'}))
            warning('ebfret:InvalidSMDFormat', ...
                    'Cannot detect SMD format from file extension. Must be one of {".mat", ".json", ".gz"}. Using ".mat".')
            args.format = 'mat';
        end
    end
    
    % initialize smd data structure
    smd = struct();
    smd.type = 'ebFRET_v_1_1_analysis';
    smd.columns = {'donor', 'acceptor', 'fret', 'viterbi_state', 'viterbi_mean'};

    % assingn global attributes
    smn.attr.num_states = analysis.dim.states;
    smd.attr.prior_mu = analysis.prior.mu;
    smd.attr.prior_beta = analysis.prior.beta;
    smd.attr.prior_W = analysis.prior.W;
    smd.attr.prior_nu = analysis.prior.nu;
    smd.attr.prior_A = analysis.prior.A;
    smd.attr.prior_pi = analysis.prior.pi;

    smd.data = struct([]);
    for n = 1:length(series)
        if ~series(n).exclude
            range = series(n).crop.min:series(n).crop.max;
            i = length(smd.data)+1;
            % store time series data
            smd.data(i).index = series(n).time(range);
            smd.data(i).values = zeros(length(range), length(smd.columns));
            smd.data(i).values(:, 1) = series(n).donor(range);
            smd.data(i).values(:, 2) = series(n).acceptor(range);
            smd.data(i).values(:, 3) = series(n).signal(range);
            smd.data(i).values(:, 4) = analysis.viterbi(n).state(range);
            smd.data(i).values(:, 5) = analysis.viterbi(n).mean(range);
            % store other analysis properties in attributes
            smd.data(i).attr.file = series(n).file;
            smd.data(i).attr.label = series(n).label;
            smd.data(i).attr.group = series(n).group;
            smd.data(i).attr.crop_min = series(n).crop.min;
            smd.data(i).attr.crop_max = series(n).crop.max;
            smd.data(i).attr.restart = analysis.restart(n);
            smd.data(i).attr.lowerbound = analysis.lowerbound(n);
            smd.data(i).attr.posterior_mu = analysis.prior.mu;
            smd.data(i).attr.posterior_beta = analysis.prior.beta;
            smd.data(i).attr.posterior_W = analysis.prior.W;
            smd.data(i).attr.posterior_nu = analysis.prior.nu;
            smd.data(i).attr.posterior_A = analysis.prior.A;
            smd.data(i).attr.posterior_pi = analysis.prior.pi;
            smd.data(i).attr.suff_stat_z = analysis.expect(n).z;
            smd.data(i).attr.suff_stat_z1 = analysis.expect(n).z1;
            smd.data(i).attr.suff_stat_zz = analysis.expect(n).zz;
            smd.data(i).attr.suff_stat_x = analysis.expect(n).x;
            smd.data(i).attr.suff_stat_xx = analysis.expect(n).xx;
        end
    end

    % run through smd.create to generate ids for each trace
    smd = ebfret.io.smd.create(...
            {smd.data.values}, ...
            smd.columns, ...
            'index', {smd.data.index}, ...
            'type', smd.type, ...
            'attr', smd.attr, ...
            'data_attrs', [smd.data.attr]);

    % save to disk
    switch args.format
        case 'mat' 
            save(args.filename, '-struct', 'smd');
        case 'json'
            ebfret.io.write_json(filename, smd, 'gzip', false);
        case 'gz'
            ebfret.io.write_json(filename, smd, 'gzip', true);
    end


