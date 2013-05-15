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

        % get prior
        u = self.analysis(a).prior;
        % populate posterior field names if necessary
        if isempty(self.analysis(a).posterior)
            self.analysis(a).posterior = u([]);
        end
        % loop over time series
        for n = args.series
            % check if we need to break
            if ~self.controls.run_analysis
                return
            end

            % get signal
            x = self.get_signal(n);

            if ~isempty(x)
                % self.set_control('series', struct('value', args.series(1)));

                % construct initial guesses for posterior parameters
                w0 = u([]);
                if (length(self.analysis(a).posterior) >= n) ...
                   && ~isempty(self.analysis(a).posterior(n).mu)
                    % always add a restart based on last result if it exists
                    w0(end+1) = self.analysis(a).posterior(n);
                end
                if args.restarts > 0
                    % first restart is uninformative
                    w0(end+1) = ebfret.analysis.hmm.init_posterior(x, u);
                end
                for r = 2:args.restarts
                    % next restarts are seeded with random draws from prior
                    w0(end+1) = ebfret.analysis.hmm.init_posterior(...
                                    x, u, ...
                                    'draw_params', true, ...
                                    'soft_kmeans', args.soft_kmeans);
                end
                
                % run variational bayes for each restart
                vb = struct();
                for r = 1:length(w0)
                    [vb(r).w vb(r).L vb(r).E] = ...
                        ebfret.analysis.hmm.vbayes(x, w0(r), u);
                end
                
                % L = arrayfun(@(vb) vb.L(end), vb);
                % [L_max r_max] = max(L);

                % determine best result
                L_max = vb(1).L(end);
                r_max = 1;
                for r = 2:length(vb)
                    if (vb(r).L(end) - L_max) > args.threshold * abs(L_max);
                        r_max = r;
                        L_max = vb(r).L(end);
                    end
                end

                % keep best result
                self.analysis(a).lowerbound(n) = vb(r_max).L(end);
                self.analysis(a).posterior(n) = vb(r_max).w;
                self.analysis(a).expect(n).z = vb(r_max).E.gamma;
                self.analysis(a).expect(n).zz = squeeze(sum(vb(r_max).E.xi, 1));
                self.analysis(a).expect(n).x = vb(r_max).E.xmean;
                self.analysis(a).expect(n).xx = vb(r_max).E.xvar + vb(r_max).E.xmean.^2;
                self.analysis(a).restart(n) = r_max + args.restarts - length(w0);

                % update plots if redraw_interval exceeded
                if toc(self.controls.redraw.series.last) > self.controls.redraw.series.interval
                    % refresh series plots
                    self.set_control('series', struct('value', n));
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
            
            % % debug:
            % L = arrayfun(@(vb) vb.L(end), vb);
            % fprintf('%s\n', sprintf('%.0e    ', (L - L_max) / abs(L_max)));

        end
        % self.refresh('ensemble');
        % drawnow();
        % last_redraw = tic();
    end
end
