function D_kl = kl_nw(w, u)
% D_kl = kl_nw(w, u)
%
% Returns Kullback-Leibler divergence 
%
%   D_kl(q || p) = Int d theta q(theta) log [q(theta) / p(theta)]
%
% for two Normal-Wishart Priors
%
%   q(mu, L | w) = Norm(mu | w.mu, w.beta L) Wish(w.W, w.nu)
%   p(mu, L | u) = Norm(mu | u.mu, u.beta L) Wish(u.W, u.nu)
%
%   p(mu, L | mu0, beta, W, nu) 
%     = B(W, nu) (beta / 2 pi)^D/2 |L|^(nu - D)/2 
%       exp[-1/2 (mu - mu0)T L (mu - mu0)]
%       exp[-1/2 Tr(Inv(W) * L)]
%
%
% Parameters
% ----------
%
%   w, u : struct 
%       Parameters for Normal-Wishart distribution
%
%       .mu (K x D)
%           Normal-Wishart prior - state means 
%       .beta (K x 1)
%           Normal-Wishart prior - state occupation count
%       .W (K x D x D)
%           Normal-Wishart prior - state precisions
%       .nu (K x 1)
%           Normal-Wishart prior - degrees of freedom
%           (must be equal to beta+1)
%   
% Output
% ------
%
%   D_kl : (K x 1)
%       Kullback-Leibler divergence 
%
%
% Jan-Willem van de Meent 
% $Revision: 1.0$ 
% $Date: 2011/08/03$

% get dimensions
[K D] = size(w.mu);

% Expectation of log emission precision |Lambda| 
% under q(W | w.W) (CB 10.65, JKC 44)
%
% E_ln_det_L(k)  =  E[ln(|Lambda(k)|)]
%                =  ln(|w.W|) + D ln(2) 
%                   + Sum_d psi((w.nu(k) + 1 - d)/2)
if D>1
    E_ln_det_L = zeros(K, 1);  
    for k=1:K
      E_ln_det_L(k) = log(det(w.W(k, :, :))) + D * log(2) + ...
                      sum(psi((w.nu(k) + 1 - (1:D)) / 2), 2);
    end
else
    E_ln_det_L = log(w.W) + D * log(2) + ...
                 sum(psi(0.5 * bsxfun(@minus, w.nu + 1, (1:D))), 2);
end

% pre-compute some terms so calculation can be vectorized
% (this was done after profiling)
if D > 1
    log_det_W_u = zeros(K, 1);
    log_det_W_w = zeros(K, 1);
    E_Tr_Winv_L = zeros(K, 1);
    dmWdm = zeros(K, 1);
    for k = 1:K
        log_det_W_u(k) = log(det(u.W(k, :, :)));
        log_det_W_w(k) = log(det(w.W(k, :, :)));
        E_Tr_Winv_L(k) = trace(inv(u.W(k, :, :) ...
                                   * w.W(k, :, :)));
        dmWdm = (w.mu(k,:) - u.mu(k,:))' ...
                * w.W(k,:,:) ...
                * (w.mu(k,:) - u.mu(k,:));
    end
else
    % for the most common D=1 case we don't need 
    % calls to det, inv, mtimes
    log_det_W_u = log(u.W);
    log_det_W_w = log(w.W);
    E_Tr_Winv_L = w.W ./ u.W;
    dmWdm = w.W .* (w.mu-u.mu).^2;
end 

% Log norm const Log[B(W, nu)] for Wishart (CB B.79)
log_B = @(log_det_W, nu) ...
        - (nu / 2) .* log_det_W ...
        - (nu * D / 2) * log(2) ...
        - (D * (D-1) / 4) * log(pi) ...
        - sum(gammaln(0.5 * bsxfun(@minus, nu + 1, (1:D))), 2);

% E_q[q(mu, L | w)]
% =
% 1/2 E_q[log |L|]
% + D log(w.beta / (2*pi)) 
% - 1/2 D
% + log(B(w.W, w.nu))
% + 1/2 (w.nu - D - 1) * E_q[log |L|]
% - 1/2 w.nu D
E_log_NW_w = 0.5 * E_ln_det_L ...        
             + 0.5 * D * log(w.beta ./ (2*pi)) ...
             - 0.5 * D ...
             + log_B(log_det_W_w, w.nu) ...
             + 0.5 * (w.nu-D-1) .* E_ln_det_L ...
             - 0.5 * w.nu * D;

% E_q[log[Norm(mu | u.mu, u.beta L)]
% =
%   1/2 D log(u.beta / 2 pi) 
%   + 1/2 E_q[log |L|]
%   - 1/2 u.beta (D / w.beta + E_q[(mu-u.mu)^T L (mu-u.mu)])  
E_log_Norm_u = 0.5 * (D * log(u.beta / (2*pi)) ...
                      + E_ln_det_L ... 
                      - D * u.beta ./ w.beta ...
                      - u.beta .* w.nu .* dmWdm);

% E_q[log[Wish(L | u.W, u.nu)]
% =
% log(B(u.W, u.nu)) 
% + 1/2 (u.nu - D - 1) E_q[log |L|]
% - 1/2 w.nu Tr[Inv(u.W) * w.W]
E_log_Wish_u = log_B(log_det_W_u, u.nu) ...
               + 0.5 * (u.nu - D - 1) .* E_ln_det_L ...
               - 0.5 * w.nu .* E_Tr_Winv_L;

% E_q[p(mu, L | u)]
% = 
% E_q[log[Norm(mu | u.mu, u.beta L)]
% + E_q[log[Wish(L | u.W, u.nu)]
E_log_NW_u = E_log_Norm_u + E_log_Wish_u;

% Dkl(q(mu, L | w) || p(mu, L | u)) 
% =
% E_q[log q(mu, L | w)] - E_q[log p(mu, L | u)]
D_kl = E_log_NW_w - E_log_NW_u;    

