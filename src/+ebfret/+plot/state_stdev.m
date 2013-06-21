function lines = state_stdev(a, b, varargin)
    % lines = state_stdev(a, b, varargin)
    % 
    % Generates plot data for the marginal on the observation noise
    % p(sigma | a, b) for a Gamma distributed prior/posterior p(l | a,b)
    % on the precision l = sigma^-2
    %
    % Inputs
    % ------
    % a, b : [K 1]
    %   Hyperparameters for Gamma distribution
    %
    % Variable Inputs
    % ---------------
    % 'xdata': cell
    %   Range of state precision values l
    % varargin : {'property', {values}}
    %   Any additional line properties.  
    %
    % Outputs
    % -------
    % lines : [K 1] struct
    %   Plot data
    %   xdata : [I 1]
    %       Range of state mean values l(i)
    %   ydata : [I 1]
    %       Probability density p(l(i) | a(k), b(k))
    %   <property> : {<values>} 
    %       Any other line properties (see doc line_props). Entries 
    %       can contain either a single or  
    %       
    ip = inputParser();
    ip.StructExpand = true;
    ip.KeepUnmatched = true;
    ip.addParamValue('labels', {}, @iscell);       
    ip.parse(varargin{:});
    args = ip.Results;

    args = ip.Results;
    K = size(a,1);

    % initialize labels if necessary
    if isempty(args.labels)
        args.labels = arrayfun(@(k) sprintf('state %d', k), ...
                        1:K, ...
                        'uniformoutput', false);
    end

    % get line properties from unmatched params
    props = cat(2, fieldnames(ip.Unmatched), struct2cell(ip.Unmatched))';
    lines = struct(props{:}, 'displayname', args.labels);
    if isscalar(lines)
        lines(1:K) = lines;
    end

    % determine x axis range if not specified
    if ~isfield(lines, 'xdata')
        for k = 1:K
            E_l = a(k,:) ./ b(k,:);
            Var_l = a(k,:) ./ b(k,:).^2;
            E_sigma = E_l.^-0.5;
            Var_sigma = 0.25 * E_l.^-3 .* Var_l;
            x = linspace(0, max(E_sigma + 4 * Var_sigma.^0.5), 101)';
            lines(k).xdata = x(2:end);
        end
    end

    % generate y values for plots
    for k = 1:K
        lines(k).ydata = ...
            bsxfun(@times, ... 
                2 * lines(k).xdata(:).^-3, ...
                exp(ebfret.analysis.dist.gamma.log_pdf(...
                    lines(k).xdata(:).^-2, ...
                    a(k,:), b(k,:))));
    end
