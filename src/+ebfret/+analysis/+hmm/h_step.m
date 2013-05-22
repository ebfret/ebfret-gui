function [u, w] = hstep(w, varargin)
% u = hstep(w, varargin)
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
%
%   w : struct  
%       Updated posterior parameters relative to new prior 
%       (assuming identical sufficient statistics)

ip = inputParser();
ip.StructExpand = true;
ip.addOptional('u', struct([]), @isstruct);       
ip.addParamValue('expect', struct([]), @isstruct);       
ip.addParamValue('max_iter', 100, @isscalar);       
ip.addParamValue('threshold', 1e-5, @isscalar);       
ip.parse(varargin{:});
args = ip.Results;

% cast posterior in array form
w_mu = cat(2, w.mu);
w_beta = cat(2, w.beta);
w_a = 0.5 * cat(2, w.nu);
w_b = 0.5 ./ cat(2, w.W);
w_A = cat(3, w.A);
w_pi = cat(2, w.pi);

u_old = args.u;
u = struct();
it = 0;
kl = [];
while true 
    % run normal-gamma updates for emission model parameters
    [u.mu, u.beta, u.a, u.b] = ...
        ebfret.analysis.dist.normgamma.h_step(...
            w_mu, w_beta, w_a, w_b, 1);
    u.nu = 2 * u.a;
    u.W = 0.5 ./ u.b;
    % add dirichlet updates for transition matrix
    u.A = ebfret.analysis.dist.dirichlet.h_step(w_A);
    % add dirichlet updates for initial state probabilities
    u.pi = ebfret.analysis.dist.dirichlet.h_step(w_pi);

    % check if we need to break
    if (it >= args.max_iter) 
        break
    end
    if (it >= 1)
        kl(it) = ebfret.analysis.hmm.kl_div(u, u_old);
    end
    if (it >= 2) & (kl(it-1) - kl(it)) ./ (1 - kl(it)) < args.threshold
        break
    end 

    % get sufficient statistics
    if (it == 0) && ~isempty(args.expect)
        E_z = cat(2, args.expect.z);
        E_z1 = cat(2, args.expect.z1);
        E_zz = cat(3, args.expect.zz);
        E_x = cat(2, args.expect.x);
        V_x = cat(2, args.expect.xx) - E_x.^2;
        % % E_zz(k,l) = w.A(k,l) - u.A(k,l)
        % E_zz = bsxfun(@minus, w_A, u_old.A);
        % E_z1 = bsxfun(@minus, w_pi, u_old.pi);
        % % E_z(k) = w.beta(k) - u.beta(k)
        % E_z = bsxfun(@minus, w_beta, u_old.beta);
        % % E_x(k) = (w.beta(k) w.mu(k) - u.beta(k) u.mu(k)) / E_z(k)
        % E_x = bsxfun(@minus, ... 
        %         w_beta .* w_mu, ...
        %         u_old.beta .* u_old.mu) ./ E_z;
        % % V_x(k) = ((1/w.W(k) - 1/u.W(k)) 
        % %            - u.beta(k) E_z(k) / (w.beta(k)) (E_z(k) - u.mu(k)).^2) 
        % %           ./ E_z
        % V_x = 2 * (bsxfun(@minus, w_b, 0.5 ./ u_old.W) ...
        %            - (bsxfun(@times, u_old.beta, E_z) ./ w_beta ...
        %               .* bsxfun(@minus, E_x, u_old.mu).^2)) ./ E_z;
    elseif isempty(args.expect)
        break
    end

    % update posteriors
    w_pi = bsxfun(@plus, E_z1, u.pi); 
    w_A = bsxfun(@plus, E_zz, u.A);
    w_beta = bsxfun(@plus, E_z, u.beta);
    w_mu = bsxfun(@plus, E_z .* E_x, u.beta .* u.mu) ./ w_beta;
    w_a = bsxfun(@plus, 0.5 * E_z, u.a);
    w_b = 0.5 * bsxfun(@plus, ...
                    2 * u.b, ...
                    E_z .* V_x ...
                    + (bsxfun(@times, E_z, u.beta) ./ w_beta ...
                        .* bsxfun(@minus, E_x, u.mu).^2));

    % increment iteration count
    it = it+1;
    u_old = u;
end

u = rmfield(u, {'a', 'b'});
u = orderfields(u, fieldnames(w));

w_nu = 2 * w_a;
w_W = 0.5 ./ w_b;
[K N] = size(w_mu);
w = struct('mu', reshape(mat2cell(w_mu, K, ones(N,1)), [N 1]), ...
           'beta', reshape(mat2cell(w_beta, K, ones(N,1)), [N 1]), ...
           'nu', reshape(mat2cell(w_nu, K, ones(N,1)), [N 1]), ...
           'W', reshape(mat2cell(w_W, K, ones(N,1)), [N 1]), ...
           'A', reshape(mat2cell(w_A, K, K, ones(N,1)), [N 1]), ...
           'pi', reshape(mat2cell(w_pi, K, ones(N,1)), [N 1]));
w = orderfields(w, fieldnames(u));