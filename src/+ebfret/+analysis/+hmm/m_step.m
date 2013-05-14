function [w ev] = m_step(u, x, g, xi)
    % w = m_step_hmm(u, x, g, xi)
    %
    % M-step: updates parameters w for q(theta | w)
    %
    % CB 10.60-10.63 and MJB 3.54 (JKC 25), 3.56 (JKC 21). 

    % get dimensions
    [K D] = size(u.mu);
    [T] = length(x);

    % Updates for mu, beta, W, nu
    [w ev] = ebfret.analysis.dist.normwish.m_step(u, x, g);

    % Update for pi
    %
    % w.pi(k) = u.pi(k) + g(1, k) 
    w.pi = u.pi + g(1, :)'; 

    % Update for A
    %
    % w.A(k, l) = u.A(k, l) + sum_t xi(t, k, l)
    w.A = u.A + squeeze(sum(xi, 1));
