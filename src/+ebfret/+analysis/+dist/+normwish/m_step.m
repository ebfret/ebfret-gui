function [w ev] = m_step(u, x, g)
    % [w ev] = m_step(u, x, g)
    %
    % M-step: updates parameters w for q(theta | w) for a
    % Normal-Wishart prior
    %
    % CB 10.60-10.63 and MJB 3.54 (JKC 25), 3.56 (JKC 21). 

    % get dimensions
    [K D] = size(u.mu);
    [T] = length(x);

    % initialize output (ensuring same fields as u)
    w = u;

    % Calculate expectation of sufficient statistics under q(z)
    %
    %   G(k) = Sum_t gamma(t, k)
    %   xmean(k,d) = Sum_t gamma(t, k)/G(k) x(t,d) 
    %   xvar(k,d,e) = Sum_t gamma(t, k)/G(k) 
    %                           (x(t,d1) - xmean(k,d))
    %                           (x(t,d2) - xmean(k,e))
    % G(k) = Sum_t gamma(t, k)
    G = sum(g, 1)';
    G = G + 1e-10;
    % g0(t,k) = g(t,k) / G(k)
    g0 = bsxfun(@rdivide, g, reshape(G, [1 K]));
    % xmean(k,d) = Sum_t g0(t, k) x(t,d) 
    xmean = sum(bsxfun(@times, ...
                       reshape(g0', [K 1 T]), ...
                       reshape(x', [1 D T])), 3);
    % dx(k,d,t) = x(t,d) - xmean(k,d) 
    dx = bsxfun(@minus, reshape(x', [1 D T]), xmean);
    % xvar(k, d, e) = Sum_t g0(t, k) dx(k, d, t) dx(k, e, t)
    xvar = squeeze(sum(bsxfun(@times, ...
                              reshape(g0', [K 1 1 T]), ... 
                              bsxfun(@times, ...
                                     reshape(dx, [K D 1 T]), ...
                                     reshape(dx, [K 1 D T]))), 4));

    % update for beta: Add counts for each state
    w.beta = u.beta + G;

    % update for mu: Weighted average of u.mu and xmean
    %
    % w.mu(k) = (u.beta(k) * u.mu(k) + G(k) * xmean(k))
    %           / (u.beta(k) + G(k))
    if D > 1
        w.mu = bsxfun(@rdivide,  ...
                      (bsxfun(@times, u.beta, u.mu) ...
                       + bsxfun(@times, G, xmean)), ...
                      w.beta);
    else
        w.mu = (u.beta .* u.mu + G .* xmean) ./ w.beta;
    end

    % update for nu: Add counts for each state
    w.nu = u.nu + G;

    % Update for W:
    %
    % Inv(w.W(k, d, e)) = Inv(u.W(k, d, e)
    %                     + G(k) * xvar(k, d, e)
    %                     + (u.beta(k) * G(k)) / (u.beta(k) + G(k))  
    %                       * (xmean(k, d) - u.mu(k, e))
    %                       * (xmean(k, d) - u.mu(k, e)) 
    
    % dx0(k,d) = xmean(k,d) - u.mu(k,d)
    dx0 = xmean - u.mu;
    % xvar0(k, d, e) = dx0(k, d) dx0(k, e)
    xvar0 = bsxfun(@times, ...
                   reshape(dx0, [K D 1]), ...
                   reshape(dx0, [K 1 D]));

    % w.W = inv(inv(u.W)(k, d1, d2) + G(k) xsigma(k, d1, d2) 
    %           + ((u.beta(k) * G(k)) / w.beta(k)) xvar0(k, d1, d2))
    if D>1
        u.W = permute(u.W, [2 3 1]);
        xvar = permute(xvar, [2 3 1]);
        xvar0 = permute(xvar0, [2 3 1]);
        for k = 1:K
            w.W(k,:,:) = inv(inv(u.W(:,:,k)) ... 
                                 + G(k) * xvar(:,:,k) ...
                                 + (u.beta(k) * G(k)) / w.beta(k) ...
                                    * xvar0(:,:,k));
        end
    else
        w.W = 1 ./ (1 ./ u.W + G .* xvar ...
                    + (u.beta .* G) ./ w.beta .* xvar0);
    end

    ev.xmean = xmean;
    ev.xvar = xvar;
