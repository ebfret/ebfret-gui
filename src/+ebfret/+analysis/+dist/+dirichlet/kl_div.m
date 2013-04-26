function D = kl_dir(alpha_p, alpha_q)
% D = kl_dir(alpha_p, alpha_q)
%
% Returns Kullback-Leibler divergence 
%
%	D_kl(p || q) = Int d theta p(theta) log [p(theta) / q(theta)] 
%
% for two Dirichlet priors
%
% 	p(pi | alpha_p) = Dir(pi | alpha_p)
% 	q(pi | alpha_q) = Dir(pi | alpha_q)
%
% Parameters
% ----------
%
%	alpha_p, alpha_q : (SZ x K)
%		Parameters of Dirichlet distributions
%
% Output
% ------
%
%	D_kl : (SZ)
%		Kullback-Leibler divergence for each set of K parameters
%		on last dimension of alpha_p and alpha_q
%
% Jan-Willem van de Meent (modified from Matthew J. Beal)
% $Revision: 1.0$ 
% $Date: 2011/08/03$

% if inputs are vectors, reshape to 1 x K
if sum(size(alpha_p)>1) <= 1
    K = length(alpha_p(:));
    alpha_p = reshape(alpha_p, [1 K]);
    alpha_q = reshape(alpha_q, [1 K]);
end

% get summation dim
d = ndims(alpha_p);

% For each set of K parameters
% 
% D_kl = log[Gamma[sum_k alpha_p(k)] 
%            / Gamma[sum_k alpha_q(k)]]
%		 + sum_k log[Gamma(alpha_p(k)) 
%                    / Gamma(alpha_q(k)) ]
%		 + sum_k (alpha_p(k) - alpha_q(k)) 
%                (psi(alpha_p(k)) - psi(Alpha_p))
D_kl = @(alpha_p, alpha_q, d) ...
       gammaln(sum(alpha_p, d)) - gammaln(sum(alpha_q,d)) ...
       - sum(gammaln(alpha_p) - gammaln(alpha_q), d) ...
       + sum((alpha_p - alpha_q) .* ...
             bsxfun(@minus, psi(alpha_p), psi(sum(alpha_p, d))), d);

% check for shared zero elements
msk = (alpha_p == 0) & (alpha_q == 0);
if any(msk(:))
    sz = size(msk);
    if length(sz) > 2
        % reshape to L x K
        L = prod(sz) / sz(d);
        K = sz(d);
        msk = reshape(msk, [L K]);
        alpha_p = reshape(alpha_p, [L K]);
        alpha_q = reshape(alpha_q, [L K]);
    else
        L = sz(1);
        K = sz(2);
    end
    D = zeros(L,1);
    for l = 1:L
        kdxs = find(~msk(l,:));
        D(l) = D_kl(alpha_p(l, kdxs), alpha_q(l, kdxs), 2);
    end
    if length(sz) > 2
        % reshape D to original dimensions
        D = reshape(D, sz(1:d-1));
    end
else
    D = D_kl(alpha_p, alpha_q, d);
end

