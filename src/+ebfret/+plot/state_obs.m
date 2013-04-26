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
    %   Mixture component index for each observation
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
    ip.addParamValue('state', [], @isvector);       
    ip.parse(varargin{:});
    args = ip.Results;

    % get assignments
    if isempty(args.state)
        args.state = ones(size(x(:)));
    end
    if isempty(args.weights)
        k_values = unique(args.state(:))';
        args.weights = 1.0 * bsxfun(@eq, args.state(:), k_values);
    end

    % get number of states 
    K = size(args.weights, 2);

    % get line properties from unmatched params
    props = cat(2, fieldnames(ip.Unmatched), struct2cell(ip.Unmatched))';
    lines = struct(props{:});
    if isscalar(lines)
        lines(1:K) = lines;
    end

    % set x-axis values if unspecified
    for k = 1:K
        if isfield(lines, 'xdata')
            [lines(k).ydata lines(k).xdata] = ...
                ebfret.plot.whist(x, lines(k).xdata, ...
                    'weights', args.weights(:,k));
        else
            [lines(k).ydata lines(k).xdata]  = ...
                ebfret.plot.whist(x, 'weights', args.weights(:,k));
        end
    end
    % normalize 
    Y = sum(sum(cat(2,lines.ydata), 1), 2);
    for k = 1:K
        dx = mean(lines(k).xdata(2:end)-lines(k).xdata(1:end-1));
        lines(k).ydata = lines(k).ydata ./ (dx .* Y);
    end
