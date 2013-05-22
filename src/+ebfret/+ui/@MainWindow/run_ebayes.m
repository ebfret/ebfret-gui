function run_ebayes(self, varargin)
    ip = inputParser();
    ip.StructExpand = true;
    ip.addParamValue('analysis', [], @isnumeric);       
    ip.addParamValue('series', [], @isnumeric);       
    ip.addParamValue('threshold', 1e-4, @isnumeric);       
    ip.addParamValue('max_iter', 100, @isscalar);       
    ip.parse(varargin{:});
    args = ip.Results;

    if isempty(args.analysis)
        if self.controls.run_all
            args.analysis = ...
                self.controls.min_states:self.controls.max_states; 
        else
            args.analysis = ...
                self.controls.ensemble.value;
        end
    end

    if isempty(args.series)
        args.series = ...
            1:length(self.series);
    end

    % ensure gui state to running
    if ~self.controls.run_analysis
        % this is necessary to prevent the next step from 
        % starting self.run_ebayes a second time
        self.controls.run_analysis = true;
        % this is necessary to make sure button state 
        % shows that analysis is running
        self.set_control('run_analysis', true);
    end

    % loop over different analysis sets
    for a = args.analysis
        % check if analysis struct needs to be initialized
        if ~isfield(self.analysis(a), 'prior') || isempty(self.analysis(a).prior) 
            self.init_analysis(a);
        end

        % switch gui to current analysis
        self.set_control('ensemble', struct('value', a));
        self.controls.redraw.ensemble.last = tic();
        drawnow();
        
        % start empirical bayes iterations
        it = 1;
        while true
            % determine number of restarts to be used in vbayes analysis
            if it == 1
                num_restarts = self.controls.init_restarts;
            else
                num_restarts = self.controls.all_restarts;
            end

            eb_interval = round(0.1 * length(self.series));
            eb_last = 0;
            for n = 1:length(self.series)
                % do variational bayes updates for posterior
                self.run_vbayes(...
                    'analysis', a, ...
                    'series', n, ...
                    'restarts', num_restarts, ...
                    'threshold', self.controls.run_precision, ...
                    'max_iter', args.max_iter);
                % % do empirical bayes updates for prior
                % if (n - eb_last) > eb_interval
                %     w = self.analysis(a).posterior(find(~[self.series.exclude]));
                %     self.analysis(a).prior = ebfret.analysis.hmm.h_step(w);
                %     eb_last = n;
                % end
            end

            % check if max iterations reached
            if (it > args.max_iter)
               break
            end
            
            % check convergence
            L(it) = sum(self.analysis(a).lowerbound);
            if (it > 1) ...
               && ((L(it) - L(it-1)) < self.controls.run_precision * abs(L(it)))  ...
                break
            end

            if it == 1
                fprintf('it %02d   L %.5e\n', it, L(it))
            else
                fprintf('it %02d   L %.5e    dL %.2e\n', it, L(it), (L(it)-L(it-1)) / abs(L(it)));
            end

            % check if we need to break
            if ~self.controls.run_analysis
                return
            end

            % run iterative empirical bayes update
            % (assuming constant posterior statistics)
            ns = find(~[self.series.exclude]);
            u = self.analysis(a).prior;
            w = self.analysis(a).posterior(find(~[self.series.exclude]));
            E = self.analysis(a).expect(find(~[self.series.exclude]));
            self.analysis(a).prior = ebfret.analysis.hmm.h_step(w, u, 'expect', E);

            % u = self.analysis(a).prior;
            % w = self.analysis(a).posterior(find(~[self.series.exclude]));
            % self.analysis(a).prior = ebfret.analysis.hmm.h_step(w);

            % increment iteration counter
            it = it + 1;
        end
    end

    self.set_control('run_analysis', 0);
    self.refresh('ensemble', 'series');
end
