function u = hstep(w)
% u = hstep(w)
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
%       .mu (K x 1)
%           Normal-Wishart prior - state means 
%       .beta (K x 1)
%           Normal-Wishart prior - state occupation count
%       .W (K x 1 x 1)
%           Normal-Wishart prior - state precisions
%       .nu (K x 1)
%           Normal-Wishart prior - degrees of freedom
%
%
% Outputs
% -------
%
%   u : struct  
%       Hyperparameters for the prior distribution p(theta | u)
%       (same fields as w)

% run normal-gamma updates for emission model parameters
[u_mu, u_beta, u_a, u_b] = ...
    ebfret.analysis.dist.normgamma.h_step(...
        cat(2, w.mu), cat(2, w.beta), 0.5 * cat(2, w.nu), 0.5 ./ cat(2, w.W));
u.mu = u_mu;
u.beta = u_beta;
u.nu = 2 * u_a;
u.W = 0.5 ./ u_b;
% add dirichlet updates for transition matrix
u.A = ebfret.analysis.dist.dirichlet.h_step({w.A});
% add dirichlet updates for initial state probabilities
u.pi = ebfret.analysis.dist.dirichlet.h_step({w.pi});
% ensure field names in correct order
u = orderfields(u, fieldnames(w));
