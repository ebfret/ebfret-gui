function lines = time_series(x, varargin)
    % lines = time_series(x, varargin)
    % 
    % Generates a time series plot.
    %
    % Inputs
    % ------
    % x : [N 1]
    %   Observations
    % t : [N 1]
    %   Time
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
    % 'crop' : struct('min', int, 'max', int)
    %   crop range for signal
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
    ip.addOptional('t', [], @isnumeric);       
    ip.addParamValue('weights', [], @isnumeric);       
    ip.addParamValue('state', [], @isnumeric);       
    ip.addParamValue('mean', [], @isnumeric);       
    ip.addParamValue('num_states', [], @isnumeric);       
    ip.addParamValue('crop', struct([]), @isstruct);       
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

    % get crop range if specified
    if ~isempty(args.crop)        
        in = args.crop.min:args.crop.max;
        out = {1:args.crop.min args.crop.max:length(x)};
    else
        in = 1:length(x);
        out = {[] []};
    end
    
    % get line properties from unmatched params
    props = cat(2, fieldnames(ip.Unmatched), struct2cell(ip.Unmatched))';
    lines = struct(props{:});

    % set time axis if not specified
    if isempty(args.t)
        args.t = 0:(length(x)-1);
    end

    % plot included area
    if ~isempty(in)
        % plot time series
        lines(1).ydata = x(in);
        lines(1).xdata = args.t(in);
        lines(1).displayname = 'signal';

        % get mean from weights if not specified
        if isempty(args.mean) & ~isempty(args.weights)
            E_x = sum(bsxfun(@times, x(in), args.weights),1) ...
                  ./ sum(args.weights, 1);
            args.mean = E_x(args.state(:));
        end
        
        if ~isempty(args.mean)
            % plot viterbi path (if supplied)
            lines(2).ydata = args.mean;
            lines(2).xdata = lines(1).xdata;
            lines(2).displayname = 'state mean';
            % plot state markers
            for k = 1:K
                lines(k+2).xdata = lines(2).xdata(args.state==k);
                lines(k+2).ydata = lines(2).ydata(args.state==k);
                lines(k+2).marker = 'o';
                lines(k+2).linestyle = 'none';
                lines(k+2).displayname = sprintf('state %d', k);
                if isfield(lines, 'color')
                    lines(k+2).markerfacecolor = lines(k+2).color;
                end
            end
        end
    end

    % plot excluded areas
    for o = 1:length(out)
        if length(out{o}) > 1
            l = length(lines)+1;
            lines(l).ydata = x(out{o});
            lines(l).xdata = args.t(out{o});
            lines(l).linestyle = '--';
            if isfield(lines, 'color')
                lines(l).color = lines(1).color;
            end
        end
    end