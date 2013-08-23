function alpha = hstep_dir(w_alpha, varargin)
% alpha = hstep_dir(w_alpha)
%
% Empirical bayes step for Dirichlet prior. Solves equations
% 
%   E[ln theta(k)] = psi(sum(alpha)) - psi(alpha(k))
%                  = sum_n weight(n)
%                    psi(sum(w_alpha(n))) - psi(w_alpha(n))
%
% Inputs
% ------
%
%   w_alpha : (N x 1) cell array or (K x N) or (L x K x N)
%       Posterior paremeters for Dirichlet distribution.
%       Each w_alpha{n} should contain either a single vector
%       of size (1 x K) or a matrix of size (L x K)
%
% Outputs
% -------
%
%  alpha : (K x 1) or (L x K) 
%       Solved Dirichlet parameters
ip = inputParser();
ip.StructExpand = true;
ip.addOptional('alpha0', [], @isnumeric);       
ip.addParamValue('max_iter', 1000, @isscalar);       
ip.addParamValue('threshold', 1e-6, @isscalar);       
ip.addParamValue('eps', 1e-12, @isscalar);       
ip.parse(varargin{:});
args = ip.Results;

% get dimensions from w_alpha
transposed = false;
if iscell(w_alpha)
    [L K] = size(w_alpha{1});
    w_alpha = cat(3, w_alpha{:});
    N = size(w_alpha, 3);
    if (K == 1) & (L > 1)
        % assume that vector valued input means L=1, 
        % and swap K with L
        [L K] = deal(K, L);
        w_alpha = reshape(w_alpha, [L K N]);
        transposed = true;
    end
else
    if ndims(w_alpha) == 3
        [L K N] = size(w_alpha);
    elseif ndims(w_alpha) == 2
        if (nargin > 2) and length(weights(:)) == 1
            % if N == 1 assuma w_alpha has size [L K]
            [L K] = size(w_alpha)
        else    
            % in absence of weights, assume N is last dimension
            L = 1;
            [K N] = size(w_alpha);
            % make sure alpha [1 K] (we'll transpose back later)
            w_alpha = reshape(w_alpha, [1 K N]);
            transposed = true;
        end
    else
        error('ebfret.analysis.dist.Dirichlet.h_step:invalidArgs', ...
              'w_alpha may be have size [K N] or [L K N]. input has size [%s]', ...
               sprintf('%d ', size(w_alpha)));
    end
end

% initialize inital guess if not specified
if isempty(args.alpha0)
    args.alpha0 = ones(L, K);
else
    args.alpha0 = reshape(args.alpha0, [L K]);
end

% theta(l,:) ~ Dirichlet
%
% E[log theta(l,k)] 
%   = mean_n (psi(w_alpha(l,k,n)) - psi(sum_k w_alpha(l,k,n)))
E_log_q = mean(bsxfun(@minus, ...
                      psi(w_alpha + args.eps), ...
                      psi(sum(w_alpha + args.eps, 2))), 3);

% do newton iterations
it = 0;
alpha = args.alpha0;
while true
    % sum of counts
    Alpha = sum(alpha, 2);
    % gradient g(k)
    g = bsxfun(@minus, psi(Alpha), psi(alpha) - E_log_q);
    % hessian H(k,l) = z + q(k) d(k,l) 
    z = psi(1, Alpha);
    q = -psi(1, alpha);
    % dalpha = (H^-1 g)(k) = (g(k) - b(k)) / q(k)
    b = sum(g ./ q, 2) ./ (1./z + sum(1./q, 2));
    dalpha = bsxfun(@minus, g, b) ./ q;
    % set constraint: alpha + dalpha >= 1e-3 alpha
    delta = min(min((1 - 1e-3) * alpha ./ dalpha, 1) .* (dalpha > 0) + (dalpha <= 0), [], 2);

    delta(find(~sum(dalpha,2)), :) = 0;
    % break if converged
    if (max(abs(dalpha(:)) ./ (alpha(:)+eps)) < args.threshold)
        break
    end
    % break if maximum iterations reached
    if it >= args.max_iter
        warning('ebfret.analysis.dist.Dirichlet.h_step:NotConverged', ...
                'Newton solver for hyperparaters did not converge in %d iterations.', ...
                max_iter)    
        break
    else 
        alpha = alpha - bsxfun(@times, delta, dalpha);
        it = it + 1;
    end
end

% ensure output has same shape as input
if transposed 
    alpha = reshape(alpha, [K 1]);
end
