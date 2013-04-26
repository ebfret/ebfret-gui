function alpha = hstep_dir(w_alpha, weights)
% alpha = hstep_dir(w_alpha, weights)
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
%       Posetrior paremeters for Dirichlet distribution.
%       Each w_alpha{n} should contain either a single vector
%       of size (1 x K) or a matrix of size (L x K)
%
%   weights : (N x 1), optional
%       Weight for each set of posterior parameters in w_alpha
%
% Outputs
% -------
%
%  alpha : (K x 1) or (L x K) 
%       Solved Dirichlet parameters


% get dimensions from w_alpha
transposed = false;
if iscell(w_alpha)
    N = length(w_alpha);
    [L K] = size(w_alpha{1});
    w_alpha = cat(3, w_alpha{:});
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
        error('hstep_dir:invalidArgs', ...
              'w_alpha may be have size [K N] or [L K N]. input has size [%s]', ...
               sprintf('%d ', size(w_alpha)));
    end
end

% set small number > eps
EPS = 10 * eps;

% set optimization settings
threshold = 1e-6;
opts = optimset('display', 'off', 'tolX', threshold, 'tolFun', eps);

% initialize dummy weights if not specified
if nargin < 2
    weights = ones([1 1 N]) / N;
elseif length(weights(:)) ~= N
    error('hstep_dir:invalidArgs', ...
          'input weights has incorrect number of elements');
else
    % ensure weights are normalized and aligned
    weights = weights / sum(weights(:));
    weights = reshape(weights, [1 1 N]);
end

% theta(l,:) ~ Dirichlet
%
% E[log theta(l,k)] 
%   = sum_n weights(n) 
%           * (psi(w_alpha(l,k,n)) - psi(sum_k w_alpha(l,k,n)))
E_log_q = sum(bsxfun(@times, weights, ...
                         bsxfun(@minus, ...
                                psi(w_alpha + EPS), ...
                                psi(sum(w_alpha + EPS, 2)))), 3);
E_q = exp(E_log_q);

% alpha(l,:) solve system of equations
alpha = mean(w_alpha, 3);
for l = 1:L
    % ignore elements with a zero mean
    kdxs = find(alpha(l,:));
    alpha_l = alpha(l, kdxs);
    E_log_q_l = E_log_q(l, kdxs);
    E_q_l = E_q(l, kdxs);
    
    % get amplitude in right ballpark first
    root_fun = @(A) (E_log_q_l - (psi(A * alpha_l) ...
                     - psi(sum(A * alpha_l)))) .* E_q_l;
    A = lsqnonlin(root_fun, 1, 0, Inf, opts);
    alpha_l = A * alpha_l;

    % now do all components
    al_old = eps * ones(size(alpha_l));
    while kl_dir(alpha_l, al_old) > threshold
        al_old = alpha_l;
        for k = 1:length(kdxs)
            a0 = sum(alpha_l .* (k ~= 1:length(kdxs)));
            root_fun = @(alk) (E_log_q_l(k) - (psi(alk) - psi(alk + a0))); 
            alpha_l(k) = lsqnonlin(root_fun, alpha_l(k), 0, Inf, opts);
        end
    end
    
    % substitute non-zero elements back in
    alpha(l,kdxs) = alpha_l;
end

% ensure output has same shape as input
if transposed 
    alpha = reshape(alpha, [K 1]);
end
