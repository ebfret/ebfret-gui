function run_vbayes(self, varargin)
    ip = inputParser();
    ip.StructExpand = true;
    ip.addParamValue('analysis', [], @isnumeric);       
    ip.addParamValue('series', [], @isnumeric);       
    ip.addParamValue('restarts', 0, @isscalar);       
    ip.addParamValue('threshold', 1e-6, @isscalar);       
    ip.addParamValue('max_iter', 100, @isscalar);       
    ip.addParamValue('soft_kmeans', false, @isscalar);       
    ip.parse(varargin{:});
    args = ip.Results;

    if isempty(args.analysis) 
        args.analysis = self.controls.ensemble.value;
    end
    if isempty(args.series)
        args.series = 1:length(self.series);
    end
    
    % loop over different analysis sets
    for a = args.analysis
        % % switch gui to current analysis
        % self.set_control('ensemble', struct('value', a));

        % % drawnow();
        % last_refresh.ensemble = tic();
        % last_refresh.series = tic();

        % break up series into batches
        b_size = 24;
        batches = {};
        b_start = 1:b_size:length(args.series);
        for b = 1:(length(b_start)-1)
            batches{b} = b_start(b):(b_start(b+1)-1);
        end
        batches{end+1} = b_start(end):length(args.series);
        % get prior
        u = self.analysis(a).prior;
        % populate posterior field names if necessary
        if isempty(self.analysis(a).posterior)
            self.analysis(a).posterior = u([]);
        end
        % loop over time series
        for b = 1:length(batches)
            % check if we need to break
            if ~self.controls.run_analysis
                return
            end
            % environment for parfor loop
            ns = batches{b};
            x = self.get_signal(ns);
            if ~iscell(x)
                % ok this is why you don't return a variable type
                % (but code may rely upon this elsewhere)
                x = {x};
                if ns > 1
                    x{ns} = x{1};
                    x{1} = [];
                end
            end
            [w(ns)] = self.analysis(a).posterior(ns);
            % output variables for parfor loop
            posterior = {};
            expect = {};
            viterbi = {};
            lowerbound = {};
            restart = {};
            % parallel loop over batch of traces
            parfor n = ns
                if ~isempty(x{n})
                    % construct initial guesses for posterior parameters
                    w0 = u([]);
                    if ~isempty(w(n).mu)
                        % always add a restart based on last result if it exists
                        w0(end+1) = w(n);
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

                    if ebfret.analysis.hmm.valid_prior(vb(r_max).w)
                        % keep best result
                        lowerbound{n} = vb(r_max).L(end);
                        posterior{n} = vb(r_max).w;
                        restart{n} = r_max + args.restarts - length(w0);
                        % calculate viterbi path
                        [viterbi{n}.state, viterbi{n}.mean] = ...
                            ebfret.analysis.hmm.viterbi_vb(vb(r_max).w, x{n});
                        % remap expected statistics
                        expect{n}.z = sum(vb(r_max).E.gamma(2:end,:), 1)';
                        expect{n}.z1 = vb(r_max).E.gamma(1, :)';
                        expect{n}.zz = squeeze(sum(vb(r_max).E.xi, 1));
                        expect{n}.x = vb(r_max).E.xmean(:);
                        expect{n}.xx = vb(r_max).E.xvar + vb(r_max).E.xmean.^2;
                    else
                        warning('ebFRET:InvalidPosterior', ...
                                'Unable to obtain a valid VB estimate for time series %d. Resetting posterior.', n);
                        lowerbound{n} = 0;
                        posterior{n} = u;
                        expect{n} = struct('z', [], 'z1', [], 'zz', [], 'x', [], 'xx', []);
                        viterbi{n} = struct('state', [], 'mean', []);
                        restart{n} = -1;
                    end
                else
                    posterior{n} = [];
                    expect{n} = [];
                    viterbi{n} = [];
                    lowerbound{n} = [];
                    restart{n} = [];
                end
            end

            % store batch results
            ns = ns(~cellfun(@isempty, posterior(ns)));
            self.analysis(a).posterior(ns) = [posterior{:}];
            self.analysis(a).expect(ns) = [expect{:}];
            self.analysis(a).viterbi(ns) = [viterbi{:}];
            self.analysis(a).lowerbound(ns) = [lowerbound{:}];
            self.analysis(a).restart(ns) = [restart{:}];

            % update plots if redraw_interval exceeded
            if toc(self.controls.redraw.series.last) > self.controls.redraw.series.interval
                % refresh series plots
                self.set_control('series', struct('value', ns(end)));
                drawnow();
                % set last refresh
                self.controls.redraw.series.last = tic();
            end
            if toc(self.controls.redraw.ensemble.last) > self.controls.redraw.ensemble.interval
                % refresh ensemble plots
                self.set_control('ensemble', struct('value', a));
                self.refresh('ensemble');
                drawnow();
                % set last refresh
                self.controls.redraw.ensemble.last = tic();
            end
        end
        % self.refresh('ensemble');
        % drawnow();
        % last_redraw = tic();
    end
end
