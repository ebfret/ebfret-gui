function lines = state_obs(x, varargin)
    % lines = state_obs(x, varargin)
    % 
    % Calculates histogram of observations for each state
    %
    % Inputs
    % ------
    % x : [N 1]
    %   Observations
    %
    % Variable Inputs
    % ---------------
    % 'xdata' : cell
    %   Histogram bin positions, or number of bins
    % 'weights' : [T K]
    %   Mixture weights for each observation
    % 'state' : [T 1]
    %   State assignment for each observation
    % varargin : {'property', {values}}
    %   Any additional line properties.  
    %
    % Outputs
    % -------
    % lines : [K 1] struct
    %   Plot data
    %   .xdata : [I 1]
    %       Bins positions for state k
    %   .ydata : [I 1]
    %       Histogram counts fpr state k
    %   .<property> : {<values>} 
    %       Any other line properties (see doc line_props). Entries 
    %       can contain either a single or  

    % parse variable args
    ip = inputParser();
    ip.StructExpand = true;
    ip.KeepUnmatched = true;
    ip.addParamValue('weights', [], @isnumeric);       
    ip.addParamValue('state', [], @isnumeric);       
    ip.addParamValue('labels', {}, @iscell);       
    ip.addParamValue('num_states', [], @isnumeric);       
    ip.parse(varargin{:});
    args = ip.Results;

    % initialize state labels if necessary
    if isempty(args.weights) & isempty(args.state)
        args.state = ones(length(x), 1);
        args.num_states = 1;
    end
    if isempty(args.num_states) 
        if ~isempty(args.weights)
            args.num_states = size(args.weights, 2);
        else
            args.num_states = max(args.state);
        end
    end 

    % get number of states
    K = args.num_states;

    % initialize labels if necessary
    if isempty(args.labels)
        args.labels = arrayfun(@(k) sprintf('state %d', k), ...
                        1:args.num_states, ...
                        'uniformoutput', false);
    end

    % get line properties from unmatched params
    props = cat(2, fieldnames(ip.Unmatched), struct2cell(ip.Unmatched))';
    lines = struct(props{:}, 'displayname', args.labels);
    if isscalar(lines)
        lines(1:K) = lines;
    end

    % generate histogram plots
    pargs.num_states = args.num_states;
    if ~isempty(args.weights)
        pargs.weights = args.weights;
    else
        pargs.state = args.state;
    end
    if isfield(lines, 'xdata')
        [ydata xdata] = ...
            ebfret.plot.whist(x, lines(1).xdata, pargs);
    else
        [ydata xdata] = ebfret.plot.whist(x, pargs);
    end
    for k = 1:K
        lines(k).xdata = xdata(:);
        lines(k).ydata = ydata(:,k);
    end

    % normalize 
    Y = sum(sum(cat(2,lines.ydata), 1), 2);
    for k = 1:K
        dx = mean(lines(k).xdata(2:end)-lines(k).xdata(1:end-1));
        lines(k).ydata = lines(k).ydata ./ (dx .* Y);
    end

