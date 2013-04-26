function lines = state_mean(m, beta, a, b, varargin)
    % lines = state_mean(m, beta, a, b, varargin)
    % 
    % Generates plot data for the marginal p(mu | m, beta, a, b)
    % for a Normal-Gamma prior/posterior
    %
    % p(mu, l | m, beta, a, b)
    %   = Norm(mu | m, beta l) Gamma(l | a, b)
    %
    % Inputs
    % ------
    % m, beta, a, b : [K *]
    %   Hyperparameters for Normal-Gamma distribution
    %
    % Variable Inputs
    % ---------------
    % 'xdata': cell
    %   Range of state mean values mu(i)
    % varargin : {'property', {values}}
    %   Any additional line properties.  
    %
    % Outputs
    % -------
    % lines : [K 1] struct
    %   Plot data
    %   xdata : [I 1]
    %       Range of state mean values mu(i)
    %   ydata : [I 1]
    %       Probability density p(mu(i) | m(k), beta(k), a(k), b(k))
    %   <property> : {<values>} 
    %       Any other line properties (see doc line_props). Entries 
    %       can contain either a single or  
    %       

    % ensure inputs shaped [K N]
    K = size(m, 1);
    m = reshape(m, [K prod(size(m))/K]);
    beta = reshape(beta, [K prod(size(beta))/K]);
    a = reshape(a, [K prod(size(a))/K]);
    b = reshape(b, [K prod(size(b))/K]);

    lines = struct(varargin{:});
    if isscalar(lines)
        lines(1:K) = lines;
    end
    if ~isfield(lines, 'xdata')
        for k = 1:K
            E_mu = m(k,:);
            V_mu = b(k,:) ./ (beta(k,:) .* a(k,:));
            lines(k).xdata = ...
                linspace(min(E_mu - 4 * V_mu.^0.5), ...
                         max(E_mu + 4 * V_mu.^0.5), 101)';
        end
    end
    for k = 1:K
        % marginal p(mu | w) is a student t with parameters
        %
        % m = m
        % l = beta a / b 
        % nu = 2 a
        lines(k).ydata = ...
            exp(ebfret.analysis.dist.studt.log_pdf(...
                    lines(k).xdata(:), ...
                    m(k,:), ...
                    beta(k,:) .* a(k,:) ./ b(k,:), ...
                    2 * a(k,:)));
    end
