function [u, w, expect] = hstep_remap(u0, expect, varargin)

ip = inputParser();
ip.StructExpand = true;
ip.addOptional('mapping', [], @isnumeric);       
ip.addParamValue('max_iter', 100, @isscalar);       
ip.addParamValue('threshold', 1e-4, @isscalar);       
ip.parse(varargin{:});
args = ip.Results;

%import ebfret.analysis.*;

% do an approximate remap of prior
if isempty(args.mapping)
    args.mapping = 1:size(u0.mu);
end
u0.a = 0.5 * u0.nu;
u0.b = 0.5 ./ u0.W;
for km = 1:max(args.mapping)
    k = find(args.mapping == km);
    msk(km) = ~isempty(k);
    if ~isempty(k)
        u.mu(km,1) = sum(u0.mu(k) .* u0.beta(k)) ./ sum(u0.beta(k));
        u.beta(km,1) = sum(u0.beta(k));
        u.a(km,1) = sum(u0.a(k));
        u.b(km,1) = sum(u0.b(k) .* u0.a(k)) ./ sum(u0.a(k));
        u.pi(km,1) = sum(u0.pi(k));
        for lm = 1:max(args.mapping)
            l = find(args.mapping == lm);
            if ~isempty(l)
                u.A(km, lm) = sum(sum(u0.A(k,l), 1), 2);
            end
        end
    end
end
u.mu = u.mu(msk);
u.beta = u.beta(msk);
u.a = u.a(msk);
u.b = u.b(msk);
u.nu = 2 * u.a;
u.W = 0.5 ./ u.b;

ns = arrayfun(@(e) ~isempty(e.z), expect);
for km = 1:max(args.mapping)
    k = find(args.mapping == km);
    msk(km) = ~isempty(k);
    if ~isempty(k)
        E_z(km, :) = arrayfun(@(e) sum(e.z(k)), expect(ns));
        E_z1(km, :) = arrayfun(@(e) sum(e.z1(k)), expect(ns));
        if isfield(expect, 'zT')
            E_zT(km, :) = arrayfun(@(e) sum(e.zT(k)), expect(ns));
        end
        E_x(km, :) = arrayfun(@(e) ebfret.nan_to_zero(sum(e.x(k) .* e.z(k)) ./ sum(e.z(k))), expect(ns));
        E_xx(km, :) = arrayfun(@(e) ebfret.nan_to_zero(sum(e.xx(k) .* e.z(k)) ./ sum(e.z(k))), expect(ns));
        for lm = 1:max(args.mapping)
            l = find(args.mapping == lm);
            if ~isempty(l)
                E_zz(km, lm, :) = ...
                    arrayfun(@(e) sum(sum(e.zz(k, l))), expect(ns)); 
            end
        end
    end
end
E_x = E_x(msk, :);
E_z = E_z(msk, :);
E_xx = E_xx(msk, :);
E_z1 = E_z1(msk, :);
E_zz = E_zz(msk, msk, :);
V_x = E_xx - E_x.^2;
if isfield(expect, 'zT')
    E_zT = E_zT(msk, :);
end

u_old = u;
it = 1;
kl = [];
while true 
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
    kl(it) = ebfret.analysis.hmm.kl_div(u, u_old);
    if (it >= 2) & (kl(it-1) - kl(it)) ./ (1 - kl(it)) < args.threshold
        break
    end 

    if isempty(expect)
        break
    end

    % increment iteration count
    it = it+1;
    u_old = u;
end

% do last update of posteriors
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

% reformat prior
u = orderfields(u, fieldnames(u0));
u = rmfield(u, {'a', 'b'});

% reformat posterior
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
w(ns) = w;

% reformat expecation values
fields = fieldnames(expect);
[K N] = size(E_x);
if isfield(expect, 'zT')
    expect = struct('x', reshape(mat2cell(E_x, K, ones(N,1)), [N 1]), ...
                    'z', reshape(mat2cell(E_z, K, ones(N,1)), [N 1]), ...
                    'z1', reshape(mat2cell(E_z1, K, ones(N,1)), [N 1]), ...
                    'zT', reshape(mat2cell(E_zT, K, ones(N,1)), [N 1]), ...
                    'xx', reshape(mat2cell(E_xx, K, ones(N,1)), [N 1]), ...
                    'zz', reshape(mat2cell(E_zz, K, K, ones(N,1)), [N 1]));
else
    expect = struct('x', reshape(mat2cell(E_x, K, ones(N,1)), [N 1]), ...
                    'z', reshape(mat2cell(E_z, K, ones(N,1)), [N 1]), ...
                    'z1', reshape(mat2cell(E_z1, K, ones(N,1)), [N 1]), ...
                    'xx', reshape(mat2cell(E_xx, K, ones(N,1)), [N 1]), ...
                    'zz', reshape(mat2cell(E_zz, K, K, ones(N,1)), [N 1]));
end
expect = orderfields(expect, fields);
expect(ns) = expect;
