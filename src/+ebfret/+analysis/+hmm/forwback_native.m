function  [gamma, xi, ln_Z] = forwback_native(px_z, A, pi)
% [gamma, xi, ln_Z] = forwback(px_z, A, pi)
%          
% Performs forward-backward message passing for HMMs.
% (fallback implemented in native matlab code)
% 
% Inputs
% ------
%
%   px_z : (T x K) 
%       Observation likelihood p(x(t) | z(t)=k, theta) = px_z(t, k) 
%       given the latent state at each time point.
%
%   A : (K x K) 
%       Transition probabilities 
%
%         p(z(t+1)=l | z(t)=k, theta)  =  A(k, l)
%
%   pi : (K x 1) 
%       Prior probabilities for states
%
%         p(z(1)=k | theta)  =  pi(k)
%
% Outputs
% -------
%
%   gamma : (T x K)
%       Posterior probabilities for states
%         p(z(t)=k | x(1:T))  =  gamma(t, k)
%
%   xi : (T-1 x K x K)
%       Posterior joint probabilities for states
%         p(z(t+1)=l, z(t)=k | x(1:T))  =  xi(t, k, l)
%
%   ln_Z : float
%       Log normalization constant Z = p(x(1:T) | theta)
%
% Jan-Willem van de Meent 
% $Revision: 1.2$  $Date: 2011/08/08$

% get dimensions
K = size(A, 1);
T = size(px_z, 1);

% Forward backward message passing  
%                           
% alpha(t, k) = p(x(1:t), z(t)=k) / p(x(1:t))
% beta(t, k) = p(x(t+1:T) | z(t)=k) / p(x(t+1:T) | x(1:t))
% scale(t, k) = p(x(t) | x(1:t-1))
alpha = zeros(T,K);
beta = zeros(T,K);
c = zeros(T,1);

% Forward pass (with scaling)
alpha(1,:) = pi' .* px_z(1, :);
c(1) = sum(alpha(1, :));
alpha(1,:) = alpha(1, :) ./ c(1); 
for t=2:T
  % alpha(t, k) = sum_l p(x(t) | z(t)=k) A(l, k) alpha(t-1, l) 
  alpha(t,:) = (alpha(t-1, :) * A) .* px_z(t, :); 
  % c(t) = p(x(t) | x(1:t-1)) = sum_l alpha(t, l)
  c(t) = sum(alpha(t, :));
  % normalize alpha by factor p(x(t) | p(x(1:t-1)))
  alpha(t,:) = alpha(t,:) ./ c(t); 
  % note: prod(c(1:t)) = p(x(1:t))
end
  
% Backward pass (with scaling)
beta(T,:) = ones(1, K);
for t=T-1:-1:1
  % beta(t, k) = sum_l p(x(t+1) | z(t)=l) A(k, l) beta(t+1, l) 
  beta(t, :) = (beta(t+1,:) .* px_z(t+1,:)) * A' / c(t+1);  
  % note: prod(c(t+1:T)) = p(x(t+1:T) | x(1:t))
end

% Posterior probabilities for states
%
% gamma(t, k) = p(z(t) | x(1:T))
%             = p(x(1:t), z(t)) p(x(t+1:T) | z(t)) / p(x(1:T))
%             = alpha(t) beta(t) 
gamma = alpha .* beta;

% Posterior transition joint probabilities
%
% xi(t, k, l) = p(z(t)=k, z(t+1)=l | x(1:T))
%             = alpha(t, k) A(k,l) px_z(t+1, l) beta(t+1, l) / c(t+1)
pxz_b_c = bsxfun(@times, 1./c(2:T), beta(2:T, :)) .* px_z(2:T, :);
xi = bsxfun(@times, ...
            reshape(alpha(1:T-1,:), [T-1 K 1]), ...
            bsxfun(@times, ...
                   reshape(A, [1, K, K]), ...
                   reshape(pxz_b_c, [T-1 1 K])));

% Evidence
%
% Ln Z = Ln[p(x(1:T))] = Sum_t log(c(t))
ln_Z = sum(log(c));