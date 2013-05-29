function lines = time_series(x, varargin)
    % lines = time_series(x, varargin)
    % 
    % Generates a time series plot.
    %
    % Inputs
    % ------
    % x : [N 1]
    %   Observations
    %
    % Variable Inputs
    % ---------------
    % 'xdata' : cell
    %   Time for each observation
    % 'weights' : [T K]
    %   Mixture weights for each observation
    % 'state' : [T 1]
    %   Mixture component index for each observation
    % 'num_states' : int
    %   Total number of states (useful when weights are not supplied)
    % 'mean' : [T 1]
    %   State mean associated with each observation
    % varargin : {'property', {values}}
    %   Any additional line properties.  
    %
    % Outputs
    % -------
    % lines : [K 1] struct
    %   Plot data
    %   .xdata : [T 1]
    %       Bins positions for state k
    %   .ydata : [T 1]
    %       Histogram counts for state k
    %   .<property> : {<values>} 
    %       Any other line properties (see doc line_props). Entries 
    %       can contain either a single or  

    % parse variable args
    ip = inputParser();
    ip.StructExpand = true;
    ip.KeepUnmatched = true;
    ip.addParamValue('weights', [], @isnumeric);       
    ip.addParamValue('state', [], @isnumeric);       
    ip.addParamValue('mean', [], @isnumeric);       
    ip.addParamValue('num_states', [], @isnumeric);       
    ip.parse(varargin{:});
    args = ip.Results;

    % get number of states 
    if isempty(args.num_states)
        if ~isempty(args.weights)
            K = size(args.weights, 2);
        elseif ~isempty(args.state)
            K = max(args.state);
        elseif ~isempty(args.mean)
            K = length(unique(args.mean));
        else
            K = 0;
        end
    else
        K = args.num_states;
    end
    % get state index from weights if not supplied
    if isempty(args.state) & ~isempty(args.weights)
        [void args.state] = max(args.weights, [], 2);
    end

    % get weights from state index if not supplied
    if ~isempty(args.state) & isempty(args.weights)
        args.weights = bsxfun(@eq, args.state(:), 1:K);
    end

    % get mean from weights if not specified
    if isempty(args.mean) & ~isempty(args.weights)
        E_x = sum(bsxfun(@times, x(:), args.weights),1) ...
              ./ sum(args.weights, 1);
        args.mean = E_x(args.state(:));
    end
    
    % get line properties from unmatched params
    props = cat(2, fieldnames(ip.Unmatched), struct2cell(ip.Unmatched))';
    lines = struct(props{:});
    % add plot for time series
    lines(1).ydata = x(:);
    % lines(1).linestyle = '-';
    % lines(1).marker = 'none';
    % lines(1).markerfacecolor = 'none';
    if ~isfield(lines, 'xdata')
        lines(1).xdata = 0:(length(x(:))-1);
    end

    % plot viterbi path (if supplied)
    if (K > 0)
        lines(2).ydata = args.mean;
        % lines(2).linestyle = '-';
        % lines(2).marker = 'none';
        % lines(2).markerfacecolor = 'none';
        if isempty(lines(2).xdata)
            lines(2).xdata = lines(1).xdata;
        end
    end

    % plot state markers
    for k = 1:K
        lines(k+2).xdata = lines(1).xdata(args.state==k);
        lines(k+2).ydata = lines(2).ydata(args.state==k);
        lines(k+2).marker = 'o';
        lines(k+2).linestyle = 'none';
        if isfield(lines, 'color')
            lines(k+2).markerfacecolor = lines(k+2).color;
        else
            lines(k+2).markerfacecolor = 'auto';
        end
    end
