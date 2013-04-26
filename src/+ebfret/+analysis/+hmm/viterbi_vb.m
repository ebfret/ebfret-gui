function [z_hat x_hat] = viterbi_vb(w, x)
% function [z_hat x_hat] = viterbi_vb(w, x)
%
% Determines Viterbi path on time series x using posterior parameter
% estimates w. 
%
% Rather than using MAP estimates for the parameters, this variant 
% uses the VB parameter estimates that are also employed
% in the forward-backward calculation of the VBEM algorithm.
%
% Inputs
% ------
%
%   w : struct
%       Variational parameters of approximate posterior distribution 
%       for parameters q(theta | w), for each of N traces 
%
%       .A (K x K)
%           Dirichlet prior for each row of transition matrix
%       .pi (K x 1)
%           Dirichlet prior for initial state probabilities
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
%   x : (TxD)
%       Observation sequence (i.e. FRET signal)
%
% Outputs
% -------
%
%   z_hat : (Tx1)
%       Index of most likely state at every time point
%
%   x_hat : (TxD)
%       Mean emissions level of most likely state at every time point
%
%
% TODO: untested for use with D>1 time series
%
% Jan-Willem van de Meent
% $Revision: 1.00$  $Date: 2011/11/07$

% get dimensions
[K D] = size(w.mu);

% get VB estimates
[E_ln_pi, E_ln_A, E_ln_px_z] = e_step_hmm(w, x);

% calculate viterbi paths
z_hat = viterbi(E_ln_px_z, E_ln_A, E_ln_pi);

% generate idealized trace
x_hat = w.mu(z_hat, :);
