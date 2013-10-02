function [z_hat, z_max, omega] = viterbi_native(ln_px_z, ln_A, ln_pi)
% function [z_hat x_hat] = viterbi(w, x)
%
% Determines the Viterbi path of most likely states for a time series
% with emission probabilities px_z.
%
% (Fallback implemented in native Matlab code)
%
%
% Inputs
% ------
%
%   ln_px_z : (T x K) 
%       Log emission probabilities p(x(t) | z(t)=k, theta) = px_z(t, k) 
%       given the latent state at each time point.
%
%   ln_A : (K x K) 
%       Log transition probabilities 
%         p(z(t+1)=l | z(t)=k, theta)  =  A(k, l)
%
%   ln_pi : (K x 1) 
%       Log prior probabilities for states
%         p(z(1)=k | theta)  =  pi(k)
%
% Outputs
% -------
%
%   z_hat : (Tx1)
%       Index of most likely state at every time point
%
%
% TODO: untested for use with D>1 time series
%
% Jan-Willem van de Meent
% $Revision: 1.10$  $Date: 2011/11/07$

% A lot of paths have 0 probablity. Not a problem for the calculation, but
% creates a lot of warning messages.
warning('off','MATLAB:log:logOfZero')

% get dimensions
[T K] = size(ln_px_z);

% intialize outputs 
% log probabilities of previous states
omega = zeros(T, K);
% most likely previous states 
z_max = zeros(T, K);
% arbitrary value, since there is no predecessor to t=1
z_max(1, :) = 0;

% forward pass
% Compute values for timestep 1
% omega(z1) = ln(p(z1)) + ln(p(x1 | z1))
omega(1, :) = ln_pi(:)' + ln_px_z(1, :);

% omega(zt) = ln(p(xt|zt)) + max{ ln(p(zt|zt-1)) + omega(zt-1) }
% CB 13.68
for t=2:T
    for k=1:K
        [omega(t, k), z_max(t, k)] = max(ln_A(:, k)' + omega(t-1, :));
        omega(t, k) = omega(t, k) + ln_px_z(t,k);
    end
end
    
% backward pass
z_hat = zeros(T, 1);
[L z_hat(T)] = max(omega(T,:));
for t=(T-1):-1:1
    z_hat(t) = z_max(t+1, z_hat(t+1));
end

warning('on','MATLAB:log:logOfZero')
