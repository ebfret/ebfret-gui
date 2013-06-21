function lines = state_dwell(alpha, varargin)
    % lines = state_dwell(alpha, varargin)
    % 
    % Calculates plot data for state life-time distributions
    %
    % Inputs
    % ------
    % alpha : [K K]
    %   Dirichlet parameters for transition probabilities
    %
    % Variable Inputs
    % ---------------
    % 'xdata' : cell
    %   Range of state precision values l
    % varargin : {'property', {values}}
    %   Any additional line properties.  
    %
    % Outputs
    % -------
    % lines : [K 1] struct
    %   Plot data
    %   xdata : [I 1]
    %       Range of state dwell times tau(i)
    %   ydata : [I 1]
    %       Probability density p(tau(i) | alpha(k,:))
    %   <property> : {<values>} 
    %       Any other line properties (see doc line_props). Entries 
    %       can contain either a single or  
    ip = inputParser();
    ip.StructExpand = true;
    ip.KeepUnmatched = true;
    ip.addParamValue('labels', {}, @iscell);       
    ip.parse(varargin{:});
    args = ip.Results;

    K = size(alpha, 1);

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

    % get beta distribution parameters for state occupancy
    for k = 1:K
        a(k,:) = alpha(k,k,:);
    end
    b = reshape(sum(alpha,2), [K size(a,2)]) - a;

    % determine x axis range if not specified
    if ~isfield(lines, 'xdata')
        for k = 1:K
            % % rho = exp(-1/tau) = a / (a+b)
            % % log(a) - log(a+b) = -1 / tau
            % % tau = 1 / (log(a+b) - log(a))
            % E_tau = 1 ./ (log(a+b) - log(a));
            % tau = -1 ./ log(rho)
            % d tau / d rho 
            % = 1 / (rho log(rho)^2)
            % = tau^2 / rho 
            % Var(tau) = (d tau / d rho)^2 Var(rho)
            %          = tau^2 / rho Var(rho)
            E_rho = a(k,:) ./ (a(k,:) + b(k,:));
            E_tau = -1 ./ log(E_rho);
            % V_rho = a(k) * b(k) / ((a(k) + b(k)).^2 * (a(k) + b(k) + 1))
            % V_tau = E_tau^2 / E_rho * V_rho
            lines(k).xdata = ...
                exp(linspace(log(0.01 * median(E_tau,2)), log(100 * median(E_tau,2)), 101))';
        end
    end
    
    % generate y values for plots
    for k = 1:K
        lines(k).ydata = ...
            bsxfun(@times, ...
                exp(-1 ./ lines(k).xdata(:)) ./ lines(k).xdata(:).^2, ...
                exp(ebfret.analysis.dist.beta.log_pdf(...
                        exp(-1 ./ lines(k).xdata(:)), a(k,:), b(k,:))));
    end
