function [counts, bins] = whist(x, varargin)
    % plots = whist(x, varargin)
    %
    % Creates histogram of observed measurement values.
    %
    % Arguments
    % ---------
    %   x : [T 1]
    %       Observations
    %   bins : int, or vector (optional)
    %       Bin positions
    %       
    % Variable Arguments
    % ------------------
    %   'weights' : [T K]
    %       Mixture weights for each observation
    %   'state' : [T 1]
    %       Mixture component index for each observation
    %
    % Outputs
    % -------
    %   counts : [B K]
    %       Observation counts
    %   bins : [B 1]
    %       Bin centres
    ip = inputParser();
    ip.StructExpand = true;
    ip.addOptional('bins', [], @isnumeric);       
    ip.addParamValue('weights', [], @isnumeric);       
    ip.addParamValue('state', [], @isvector);       
    ip.parse(varargin{:});
    args = ip.Results;
    % get bins
    if isempty(args.bins)
        args.bins = max(10, min(length(x(:)) / 10, 200));
    end
    if isscalar(args.bins)
        [void bins] = hist(x(:), args.bins);
    else
        bins = args.bins;
    end
    bins = bins(:);
    % get assignments
    if isempty(args.state)
        args.state = ones(size(x(:)));
    end
    if isempty(args.weights)
        k_values = unique(args.state(:))';
        args.weights = 1.0 * bsxfun(@eq, args.state(:), k_values);
    end
    % get bin index for each observation
    [void, x_bin] = min(abs(bsxfun(@minus, x(:), bins(:)')), [], 2);
    for b = 1:length(bins)
        counts(b, :) = sum(args.weights(x_bin==b,:), 1);
    end
