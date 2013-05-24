function [rep par] = report(signal, prior, expect, varargin)
    ip = inputParser();
    ip.StructExpand = true;
    ip.KeepUnmatched = true;
    ip.addRequired('x', @iscell);       
    ip.addRequired('prior', @isstruct);       
    ip.addRequired('expect', @isstruct);       
    ip.addParamValue('splits', {}, @iscell);       
    ip.addParamValue('labels', {}, @iscell);       
    ip.addParamValue('mapping', [], @isnumeric);       
    ip.parse(signal, prior, expect, varargin{:});
    args = ip.Results;

    import ebfret.analysis.*
    % if nargin < 2
    %     [fname, fpath, findex] = ...
    %     uiputfile({'*.csv', 'ebFRET analysis summary (.csv)';});
    %     filename = sprintf('%s/%s', fpath, fname);
    % end

    if isempty(args.splits)
        args.splits = {1:length(args.x)};
    end
    
    for s = 1:length(args.splits)
        % set label
        % basic statistics on number of series and time points 
        ns = args.splits{s};
        x = args.x(ns);
        T = cellfun(@length, x);
        ns = ns(T>0);
        T = T(T>0);
        if ~isempty(args.labels)
            rep(s).Series.Label = args.labels{s};
        end
        rep(s).Series.Number = length(T);
        rep(s).Series.Length.Mean = mean(T);
        rep(s).Series.Length.Std = std(T);
        rep(s).Series.Length.Median = median(T);
        rep(s).Series.Length.Min = min(T);
        rep(s).Series.Length.Max = max(T);
        rep(s).Series.Length.Total = sum(T);
        % get prior posterior and expectation values for splits
        [u, w, e] = hmm.h_step_remap(...
                        args.prior, args.expect(ns), args.mapping);
        par(s).prior = u;
        par(s).posterior(ns) = w;
        par(s).expect(ns) = e;
        % statistics
        if ~isempty(args.labels)
            rep(s).Statistics.Label = {args.labels{s}};
            for k = 2:length(u.mu)
                rep(s).Statistics.Label{k} = '';
            end
        end
        rep(s).Statistics.State = (1:length(u.mu))';
        rep(s).Statistics.Occupancy.Mean = ...
            normalize(sum(cat(2, e.z), 2));
        rep(s).Statistics.Occupancy.Total = ...
            sum(cat(2, e.z), 2);
        rep(s).Statistics.Observation.Mean = ...
            sum(cat(2, e.x) .* cat(2, e.z), 2) ./ sum(cat(2, e.z), 2);
        rep(s).Statistics.Observation.Std = ...
            (sum(cat(2, e.xx) .* cat(2, e.z), 2) ./ sum(cat(2, e.z), 2) ...
             - rep(s).Statistics.Observation.Mean.^2).^0.5;
        rep(s).Statistics.Transitions.Mean = ...
            normalize(sum(cat(3, e.zz), 3),2);
        rep(s).Statistics.Transitions.Total = ...
            sum(cat(3, e.zz), 3);

        % parameter values
        if ~isempty(args.labels)
            rep(s).Parameters.Label = {args.labels{s}};
            for k = 2:length(u.mu)
                rep(s).Parameters.Label{k} = '';
            end
        end
        rep(s).Parameters.State = (1:length(u.mu))';
        [rep(s).Parameters.Center.Mean rep(s).Parameters.Precision.Mean] = ...
            dist.normgamma.mean(u.mu, u.beta, 0.5 * u.nu, 0.5 ./ u.W);
        [V_mu, V_l] = ...
            dist.normgamma.var(u.mu, u.beta, 0.5 * u.nu, 0.5 ./ u.W);
        [rep(s).Parameters.Center.Std rep(s).Parameters.Precision.Std] = ...
            deal(V_mu.^0.5, V_l.^0.5);
        [void rep(s).Parameters.Precision.Mode] = ...
            dist.normgamma.mode(u.mu, u.beta, 0.5 * u.nu, 0.5 ./ u.W);
        rep(s).Parameters.Dwell_Time.Mode = ...
            dist.dirichlet.tau(u.A);
        % rep(s).variance.mean = ...
        %     dist.invgamma.mean(0.5 * u.nu, 0.5 * u.W);
        % rep(s).variance.mode = ...
        %     dist.invgamma.mode(0.5 * u.nu, 0.5 * u.W);
        % rep(s).variance.var = ...
        %     dist.invgamma.var(0.5 * u.nu, 0.5 * u.W);
        rep(s).Parameters.Transition_Matrix.Mean = ...
            dist.dirichlet.mean(u.A);
        % rep(s).Parameters.Transition_Matrix.Mode = ...
        %     dist.dirichlet.mode(u.A);
        rep(s).Parameters.Transition_Matrix.Std = ...
            dist.dirichlet.var(u.A).^0.5;
    end
end