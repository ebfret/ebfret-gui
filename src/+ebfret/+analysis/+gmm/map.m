function [theta, L] = map(x, theta0, u, varargin)
% [theta, L] = map(x, theta0, varargin)
%
% Maximum Posterior Expectation Maximization for a 
% Gaussian Mixture Model. 
%
%
% Inputs
% ------
%
%   x : (T x D) float
%     Observation time series. May be one dimensional 
%     (e.g. a FRET signal), or higher dimensional 
%     e.g. D=2 for a donor/acceptor signal)
%
%   theta0 : struct 
%     Initial guess for the model parameters.
%
%     .pi (K x 1)
%         Initial state probabilities
%     .mu (K x D)
%         State emissions means 
%     .Lambda (K x D)
%         State precisions
%
%   u : struct
%     Prior parameters
%
%     .pi : (K x 1)
%         Dirichlet prior for initial state probabilities
%     .mu : (K x D)
%         Normal-Wishart prior - state means 
%     .beta : (K x 1)
%         Normal-Wishart prior - state occupation count
%     .W : (K x D x D)
%         Normal-Wishart prior - state precisions
%     .nu : (K x 1)
%         Normal-Wishart prior - degrees of freedom
%         (should be equal to beta+1)
%
%
% Variable Inputs
% ---------------
%
%   threshold : float (default: 1e-5)
%      Convergence threshold. Execution halts when the relative
%      increase in the lower bound evidence drop below threshold 
%
%   max_iter : int (default: 100)
%      Maximum number of iteration before execution is truncated
%
%
% Outputs
% -------
%
%   theta : struct
%       Maximum likelihood value for parameters (same fields as theta0)
%
%   L : (I x 1) float
%       Log likelihood for each iteration
%
%
% Jan-Willem van de Meent
% $Revision: 1.00 $  $Date: 2012/06/09$

Debug = false;

% parse inputs
ip = inputParser();
ip.StructExpand = true;
ip.addRequired('x', @(x) isnumeric(x) & (ndims(x)==2));
ip.addRequired('theta0', @(u) all(isfield(u, {'pi', 'mu', 'Lambda'})));
ip.addRequired('u', @(u) all(isfield(u, {'pi', 'mu', 'beta', 'W', 'nu'})));
ip.addParamValue('threshold', 1e-5, @isscalar);
ip.addParamValue('max_iter', 100, @isscalar);
ip.parse(x, theta0, u, varargin{:});

% collect inputs
args = ip.Results;
x = args.x;
theta0 = args.theta0;

% set theta to initial guess
theta = theta0;

% get dimensions
[T D] = size(x);
K = length(theta.pi);

% precompute log determinant and inverse of W (need these later)
for k = 1:K
    W = squeeze(u.W(k,:,:));
    log_det_W(k) = log(det(W));
    W_inv(:, :, k) = inv(W);
end

% log norm const log(C(alpha)) for Dirichlet (CB B.23)
log_C = gammaln(sum(u.pi)) - sum(gammaln(u.pi));

% log norm const log[B(W, nu)] for Wishart (CB B.79)
log_B = - (u.nu / 2) .* log_det_W(:) ...
        - (u.nu * D / 2) * log(2) ...
        - (D * (D-1) / 4) * log(pi) ...
        - sum(gammaln(0.5 * bsxfun(@minus, u.nu + 1, (1:D))), 2);

% Main loop of algorithm
for it = 1:args.max_iter
    % E-STEP: UPDATE Q(Z)
    if D>1
        % Calculate Mahalanobis distance
        Lambda = permute(theta.Lambda, [2 3 1]);
        % dx(d, t, k) = x(t, d) - mu(k, d)
        dx = bsxfun(@minus, x', reshape(theta.mu', [D 1 K]));
        % dxLdx(t, k) = Sum_de dx(d,t,k) * l(d, e, k) * dx(e, t, k)
        Delta2 = squeeze(mtimesx(reshape(dx, [1 D T K]), ...
                                 mtimesx(reshape(Lambda, [D D 1 K]), ...
                                         reshape(dx, [D 1 T K]))));
        % note, the mtimesx function applies matrix multiplication to
        % the first two dimensions of an N-dim array, while using singleton
        % expansion to the remaining dimensions. 
        % 
        % TODO: make mtimesx usage optional? (needs compile on Linux/MacOS)

        % calculate precision matrix determinant
        ln_det_Lambda = zeros([K 1]);
        for k = 1:K
            ln_det_Lambda(k) = log(det(Lambda(:,:,k)));
        end
    else
        % dx(t, k) = x(t) - mu(k)
        dx = bsxfun(@minus, x, theta.mu');
        % dxLdx(t, k) = Sum_de dx(t,k) * L(k) * dx(t, k)
        Delta2 = bsxfun(@times, dx, bsxfun(@times, theta.Lambda', dx));
        % get precision 'determinant'
        ln_det_Lambda = log(theta.Lambda);
    end

    % Log emissions probability p(x | z, theta) 
    %
    % ln_px_z(t, k)
    %   = log(1 / 2 pi) * (D / 2)
    %     + 0.5 * E[ln(|Lambda(k,:,:)|]]
    %     - 0.5 * E[Delta(t,k)^2]
    ln_px_z = log(2 * pi) * (-D / 2) ...
                + bsxfun(@plus, 0.5 * ln_det_Lambda', ...
                               -0.5 * Delta2);

    % subtract mean from ln_px_z (avoids under/overflow in next step)
    ln_px_z0 = bsxfun(@minus, ln_px_z, mean(ln_px_z, 2));

    % Calculate responsibilities g = q(z) = p(z | x, theta)
    %
    % g(t, k) = px_z(t, k) * pi(k) 
    %           / sum_l px_z(t, l) * theta.pi(l)
    g = ebfret.normalize(bsxfun(@times, exp(ln_px_z0), theta.pi(:)'), 2);

    % Calculate log likelihood
    %
    % L_ml = sum_{t,k} g(t, k) log(px_z(t, k) * pi(k))
    ln_pxz = bsxfun(@plus, ln_px_z, log(theta.pi(:)'));
    L_ml = sum(g(:) .* ln_pxz(:));

    % Calculate log prior
    %
    % log(p(theta)) = log(p(pi)) + log(p(mu)) + log(p(Lambda))

    % log(p(pi | u.pi)) (Dirichlet)
    ln_pq_pi = log_C + sum(theta.pi .* (u.pi - 1));

    % log(p(mu | u.mu, u.beta * Lambda)) (Multivariate Normal)
    if D == 1
        Delta2 = u.beta(:) .* theta.Lambda(:) .* (theta.mu(:) - u.mu(:)).^2;
    else
        dmu = theta.mu - u.mu;
        for k = 1:K
            Delta2(k) = dmu(k, :) * (u.beta(:) .* squeeze(theta.Lambda(k,:,:))) * dmu(k,:)';
        end
    end
    ln_pq_mu = log(2 * pi) * (-D / 2) ...
               + 0.5 * (u.beta(:) * D + ln_det_Lambda(:)) ...
               - 0.5 * Delta2(:);

    % log(p(Lambda | u.W, u.nu)) (Wishart)
    if D == 1 
        Tr_Winv_L = W_inv(:) .* theta.Lambda(:);
    else
        for k = 1:K
            Tr_Winv_L(k) = trace(W_inv(:,:,k) * squeeze(theta.Lambda(:,:,k)));
        end
    end
    ln_pq_L = log_B(:) + ...
              0.5 * ln_det_Lambda(:) .* (u.nu(:) - D - 1) ...
              - 0.5 * Tr_Winv_L(:);

    % Calculate log joint log(p(x, theta))
    %
    % L = L_ml + log(p(theta | u))
    L(it) = L_ml + ln_pq_pi + sum(ln_pq_mu) + sum(ln_pq_L);

    % issue warning if log joint decreases
    if it>2 && ((L(it) - L(it-1)) < -10 * args.threshold * abs(L(it))) 
        warning('gmm_map:joint_decreased', ...
                'Log joint p(x,theta) decreased by %e', ...
                L(it) - L(it-1));
        % fprintf('L_ml: %.2e, ln_pq_pi: %.2e, ln_pq_mu: %.2e, ln_pq_L: %.2e\n', ...
        %         L_ml, ln_pq_pi, sum(ln_pq_mu), sum(ln_pq_L));
    end

    % check if the lower bound increase/decrease is less than threshold
    if (it>2)    
        if abs((L(it) - L(it-1)) / L(it-1)) < args.threshold || ~isfinite(L(it)) 
            L(it+1:end) = [];  
            break;
        end
    end

    % M STEP: UPDATE THETA

    % calculate posterior parameters (same  as VBEM m step)
    w = ebfret.analysis.dist.normwish.m_step(u, x, g);

    % get map values
    theta.pi = ebfret.normalize(w.beta);
    theta.mu = w.mu;
    theta.Lambda = bsxfun(@times, w.W, w.nu);
end