function [posterior, expect, lowerbound, restart] = vbayes(x, prior, varargin)
    % analysis = ebayes(series, varargin)
    %
    % Run empirical bayes analysis batch job
    ip = inputParser();
    ip.StructExpand = true;
    ip.addParamValue('posterior', []);       
    ip.addParamValue('restarts', 0, @isscalar);       
    ip.addParamValue('threshold', 1e-6, @isscalar);       
    ip.addParamValue('max_iter', 100, @isscalar);       
    ip.addParamValue('soft_kmeans', false, @isscalar);       
    ip.parse(varargin{:});
    args = ip.Results;

    % output variables for parfor loop
    posterior = {};
    expect = {};
    viterbi = {};
    lowerbound = {};
    restart = {};

    % parallel loop over traces
    u = prior;
    parfor n = 1:length(x)
        if ~isempty(x{n})
            % construct initial guesses for posterior parameters
            w0 = u([]);
            try 
                % always add a restart based on last result if it exists
                if ~isempty(args.posterior(n).mu)
                    w0(end+1) = args.posterior(n);
                end
            catch
                % if posterior not supplied this fails (but no biggie)
            end
            if args.restarts > 0
                % first restart is uninformative
                w0(end+1) = ebfret.analysis.hmm.init_posterior(x{n}, u);
            end
            for r = 2:args.restarts
                % next restarts are seeded with random draws from prior
                w0(end+1) = ebfret.analysis.hmm.init_posterior(...
                                x{n}, u, ...
                                'draw_params', true, ...
                                'soft_kmeans', args.soft_kmeans);
            end
            
            % run variational bayes for each restart
            vb = struct();
            for r = 1:length(w0)
                [vb(r).w vb(r).L vb(r).E] = ...
                    ebfret.analysis.hmm.vbayes(x{n}, w0(r), u);
            end
            
            % determine best result
            L_max = vb(1).L(end);
            r_max = 1;
            for r = 2:length(vb)
                if (vb(r).L(end) - L_max) > 1e-2 * args.threshold * abs(L_max) ...
                    || isnan(L_max)
                    r_max = r;
                    L_max = vb(r).L(end);
                end
            end

            % keep best result
            lowerbound{n} = vb(r_max).L(end);
            posterior{n} = vb(r_max).w;
            restart{n} = r_max + args.restarts - length(w0);

            % remap expected statistics
            expect{n}.z = sum(vb(r_max).E.gamma(2:end,:),1)';
            expect{n}.z1 = vb(r_max).E.gamma(1,:)';
            expect{n}.zz = squeeze(sum(vb(r_max).E.xi, 1));
            expect{n}.x = vb(r_max).E.xmean(:);
            expect{n}.xx = vb(r_max).E.xvar + vb(r_max).E.xmean.^2;
        else
            posterior{n} = [];
            expect{n} = [];
            viterbi{n} = [];
            lowerbound{n} = [];
            restart{n} = [];
        end
    end

    % reformat output to struct array
    posterior = [posterior{:}];
    expect = [expect{:}];
    viterbi = [viterbi{:}];
    lowerbound = [lowerbound{:}];
    restart = [restart{:}];
end