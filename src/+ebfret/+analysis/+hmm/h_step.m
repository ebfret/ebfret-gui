function u = hstep_hmm(w, weights)
% u = hstep_hmm(w, weights)
%
% Hyper parameter updates for empirical Bayes inference (EB)
% on a single-molecule FRET dataset.
%
% The input of this method is a set of  posteriors produced by running 
% Variational Bayes Expectation Maxmimization (VBEM) on a series
% of FRET traces and maximizes the total summed evidence by solving 
% the system of equations:
%
%   Sum_n  Grad_u L_n  =  0
%
% Where u is the set of hyperparameters that determine the form of
% the prior distribution on the parameters. L_n is the lower bound
% evidence for the n-th trace, defined by:
%
%   L  =  Integral d z d theta  q(z) q(theta)  
%         ln[ p(x, z | theta) p(theta) / (q(z) q(theta)) ]
%
% The approximate posteriors q(z) and q(theta), which have been
% optimised in the VBEM process, are now kept constant as Sum L_n is
% maximised wrt to u.
%
% 
% Inputs
% ------
%
%   w : struct (N x 1)
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
%
%   weights : (N x 1) optional
%       Weighting for each trace in updates.
%
%
% Outputs
% -------
%
%   u : struct  
%       Hyperparameters for the prior distribution p(theta | u)
%       (same fields as w)
%
% Jan-Willem van de Meent
% $Revision: 1.10$  $Date: 2011/08/04$

% Note: this currently does not work for 2D donor/acceptor inference
%
% Explanation of Updates
% ----------------------
%
% For a prior/posterior in Conjugate Exponential form 
%
%   p(theta)  =  f(nu, chi) g(theta)^nu exp[ eta(theta) . chi ]
%   q(theta)  =  f(nu', chi') g(theta)^nu' exp[ eta(theta) . chi' ]
%
% The derivatives of L w.r.t. to the parameters are given by
% 
%   Grad_nu L  =  (Grad_nu f(eta, nu)) / f(eta, nu)
%                 + <ln g(theta)>_q(theta)  
%
% Similarly the derivatives wrt to nu_i are given by                   
%
%   Grad_chi_i L =  (Grad_chi_i f(eta, nu)) / f(eta, nu)
%                   + <eta_i>_q(theta)
%
% We now define the averaged expectation value over all traces:
%
%   E[theta]  =  1/N  Sum_n  <theta>_q(theta | w(n))
%
% The system of equations Grad Sum_n L^n = 0 can now be solved
% from the above identities, resulting in the following updates
%
%   (mu, lambda) ~ NormalGamma(mu, lambda | u.mu, u.beta, u.nu, u.W)
%
%   u.beta  =  u.nu - 1
%   u.mu  =  E[mu lambda] / E[lambda]
%   u.W  =  E[lambda] / u.nu
%   -<log g> =  1/2 ( 1 / (u.nu -1)
%                     + E[mu lambda]^2 / E[lambda]
%                     + Log[pi / uW]
%                     - psi(u.nu/2) )
%       
% {A(k,:)} ~ Dir(A_k | u.A(k,:))
%
%   psi(Sum u.A(k,:)) - psi(u.A(k,:)) = -E[log A(k,:)]
%
% pi ~ Dir(pi | u.pi)
%
%   psi(Sum u.pi) - psi(u.pi) = -E[log pi] 

% intialize empty weights if unspecified
if nargin < 2
    weights = ones(size(w));
end
% run normal-wishart updates for emission model parameters
u = ebfret.analysis.dist.normwish.h_step(w, weights);
% add dirichlet updates for transition matrix
u.A = ebfret.analysis.dist.dirichlet.h_step({w.A}, weights);
% add dirichlet updates for initial state probabilities
u.pi = ebfret.analysis.dist.dirichlet.h_step({w.pi}, weights);
% ensure field names in correct order
u = orderfields(u, fieldnames(w));
