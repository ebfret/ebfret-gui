function [u_m, u_beta, u_a, u_b] = h_step(w_m, w_beta, w_a, w_b, varargin)
% [u_m, u_beta, u_a, u_b] = h_step(w_m, w_beta, w_a, w_b, varargin)
%
% Emperical Bayes step for a Normal-Gamma prior. Solves equations
%
%     m = - E[m l] /  E[-l]
%     beta = 1 / (E[m^2 l] - E[m l]^2 / E[l])
%     psi(a) - log(a) = E[log l] - log(E[l])
%     b = - a / E[-l]
% 
% Inputs
% ------
%
% w_m, w_beta, w_a, w_b : [K N]
%     Normal-Gamma posterior parameters for K states and N samples
%
% Outputs
% -------
%
% u_m, u_beta, u_a, u_b : [K 1]
%     Estimated hyperparameters

ip = inputParser();
ip.StructExpand = true;
ip.addOptional('a0', 2, @isnumeric);       
ip.addParamValue('max_iter', 1000, @isscalar);       
ip.addParamValue('threshold', 1e-6, @isscalar);       
ip.addParamValue('eps', 1e-12, @isscalar);       
ip.addParamValue('constraints', ...
    struct('beta', [1e-3,inf], 'a', [1+1e-3, inf]), @isscalar);       
ip.parse(varargin{:});
args = ip.Results;

% accept cell arrays
if iscell(w_m)
    w_m = cat(2, w_m{:});
end
if iscell(w_beta)
    w_beta = cat(2, w_beta{:});
end
if iscell(w_a)
    w_a = cat(2, w_a{:});
end
if iscell(w_b)
    w_b = cat(2, w_b{:});
end

% initialize u from naive average
u_m = mean(w_m, 2);
u_beta = mean(w_beta, 2);
u_a = mean(w_a, 2);
u_b = mean(w_b, 2);

% Expectation Value E[lambda] = w_a / w_b)
E_l = mean(w_a ./ w_b, 2);
log_E_l = log(E_l);

% Expectation Value E[mu * lambda] = w_m * w_a / w_b
E_ml = mean(w_m .* w_a ./ w_b, 2);

% Expectation Value E[mu^2 * lambda] = 1 / w_beta + w_m^2 w_a / w_b
E_m2l = mean(1 ./ w_beta + w_m.^2 .* w_a ./ w_b, 2);

% Expectation Value E[log lambda] = psi(w_a) - log(w_b)
E_log_l = mean(psi(w_a) - log(w_b), 2);

% solve for u_a
%
%   psi(u_a) - log(u_a) 
%       = E[log(lambda)] - log(E[lambda])
a_old = args.a0;
a = args.a0;
it = 0;
while true
    % break if maximum iterations reached
    if it >= args.max_iter
        warning('ebfret.analysi.dist.normgamma.h_step:NotConverged', ...
                'Newton solver for hyperparaters did not converge in %d iterations.', ...
                args.max_iter)    
        break
    end
    % gradient 
    % g = psi(a) - log(a) - (E[log(l)] - log(E[l]))
    g = psi(a) - log(a) - (E_log_l - log_E_l);
    % hessian H 
    H = psi(1,a) - 1 ./ a; 
    % da = (H^-1 g)
    da = g ./ H;
    % ensure a > 1 
    a = max(a - da, 1 + 1e-3 * (a - 1));
    % break if converged
    if all((abs(a - a_old) ./ a) < args.threshold)
        break
    end
    it = it + 1;
    a_old = a;
end

% m = E[m l] / E[l]
u_m = E_ml ./ E_l;

% beta = 1 / (E[m^2 l] - E[m l]^2 / E[l])
u_beta = 1 ./ (E_m2l - E_ml.^2 ./ E_l);

% psi(u_a) - log(u_a) = E[log(lambda)] - log(E[lambda])
u_a = a;

% b = a / E[l]
u_b = a ./ E_l;

% enforce constraints
vars = fieldnames(args.constraints);
u = struct('m', u_m, 'beta', u_beta, 'a', u_a, 'b', u_b);
for v = 1:length(vars)
    if isfield(u, vars{v})
        msk = u.(vars{v}) < args.constraints.(vars{v})(1);
        u.(vars{v})(msk) = args.constraints.(vars{v})(1);
        msk = u.(vars{v}) > args.constraints.(vars{v})(2);
        u.(vars{v})(msk) = args.constraints.(vars{v})(2);
    end
end

[u_m, u_beta, u_a, u_b] = deal(u.m, u.beta, u.a, u.b);