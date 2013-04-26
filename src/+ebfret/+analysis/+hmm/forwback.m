function  [gamma, xi, ln_Z] = forwback(px_z, A, pi)
% [gamma, xi, ln_Z] = forwback(px_z, A, pi)
%          
% Performs forward-backward message passing for HMMs.
% (implemented in C)
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

% SEE forwback.c for implementation
