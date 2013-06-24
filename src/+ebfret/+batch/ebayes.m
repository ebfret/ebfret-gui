function [analysis, options] = ebayes(series, options, clip)
    % analysis = ebayes(series, options, clip)
    %
    % Run empirical bayes analysis job(s).

    % clip and crop time series
    if nargin < 3
        clip.min = -inf;
        clip.max = inf;
    end
    for n = 1:length(series)
        s = series(n);
        if ~s.exclude
            x{n} = s.signal(s.crop.min:s.crop.max);
            % clip to specified range
            if clip.min < inf
                tclip = find(x{n} < clip.min);
                x{n}(tclip) = clip.min;
            end
            if clip.max < inf
                tclip = find(x{n} > clip.max);
                x{n}(tclip) = clip.max;
            end
        else
            x{n} = [];
        end
    end

    for o = 1:length(options(:))
        opt = options(o);
        % initialize prior
        if ~isfield(opt, 'prior') || isempty(opt.prior.mu)
            if ~isfield(opt, 'init')
                % guess prior
                try 
                    K = opt.dim.states;
                catch 
                    K = b;
                end
                analysis(o).prior = ebfret.analysis.hmm.guess_prior(x, K);
            else
                % initialize prior from specified settings
                analysis(o).prior = ebfret.analysis.hmm.init_prior(...
                                        opt.init.theta, opt.init.counts);
            end
            analysis(o).dim.states = length(analysis(o).prior.mu);
            % store resulting prior in opt
            opt.prior = analysis(o).prior;
        else 
            % use specified initial guess for prior
            analysis(o).prior = opt.prior;
        end
        % store number of states
        K = length(analysis(o).prior.mu);
        analysis(o).dim.states = K; 
        opt.dim.states = K; 

        % initialize convergence options
        if ~isfield(opt, 'eb')
            opt.eb = struct();
        end
        if ~isfield(opt.eb, 'max_iter')
            opt.eb.max_iter = 100;
        end
        if ~isfield(opt.eb, 'threshold')
            opt.eb.threshold = 1e-4;
        end
        if ~isfield(opt, 'vb')
            opt.vb = struct();
        end
        if ~isfield(opt.vb, 'max_iter')
            opt.vb.max_iter = opt.eb.max_iter;
        end
        if ~isfield(opt.vb, 'threshold')
            opt.vb.threshold = 0.1 * opt.eb.threshold;
        end
        if ~isfield(opt.vb, 'restarts')
            opt.vb.restarts = 2;
        end

        % analysis
        it = 1;
        analysis(o).posterior = analysis(o).prior([]);
        while true
            % run variational bayes on each time series
            [posterior, expect, lowerbound] = ...
                ebfret.batch.vbayes(x, analysis(o).prior, ...
                    'posterior', analysis(o).posterior, opt.vb);
            ns = find(~[series.exclude]);
            analysis(o).posterior(ns) = posterior;
            analysis(o).expect(ns) = expect;
            analysis(o).lowerbound(ns) = lowerbound;

            % calculate and print out sum lower bound
            L(it) = sum(analysis(o).lowerbound);
            if it == 1
                fprintf('K%02d    it %02d   L %.5e\n', K, it, L(it))
            else
                fprintf('K%02d    it %02d   L %.5e    dL %.2e\n', K, it, L(it), (L(it)-L(it-1)) / abs(L(it)));
            end

            % check if max iterations reached
            if (it > opt.eb.max_iter)
               break
            end
            
            % check convergence
            if (it > 1) ...
               && ((L(it) - L(it-1)) < opt.eb.threshold * abs(L(it)))  ...
                break
            end

            % run iterative empirical bayes update
            % (assuming constant posterior statistics)
            analysis(o).prior = ebfret.analysis.hmm.h_step(...
                                    analysis(o).posterior, ...
                                    analysis(o).prior, ...
                                    'expect', analysis(o).expect);

            % increment iteration counter
            it = it + 1;
        end

        % calculate viterbi paths
        ns = find(~[series.exclude]);
        for n = ns
            [analysis(o).viterbi(n).state, analysis(o).viterbi(n).mean] = ...
                ebfret.analysis.hmm.viterbi_vb(analysis(o).posterior(n), x{n});
        end

        % store 
        opt.L = L;
        opts(o) = opt;
    end
    analysis = reshape(analysis, size(options));
    options = reshape(opts, size(options));
end