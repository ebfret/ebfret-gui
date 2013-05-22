function [M_tau] = tau(alpha)
    % M_tau = tau(alpha)
    % 
    % Most likely state dwell time
    import ebfret.analysis.dist.*
    if nargin < 2
        precision = 100*eps;
    end
    transposed = isrow(alpha);
    if transposed
        alpha = alpha';
    end
    a = diag(alpha);
    b = sum(alpha,2) - a;
    log_B = gammaln(a) + gammaln(b) - gammaln(a+b);
    E_rho = a ./ (a+b);
    for k = 1:length(a)
        % calculate expectation and variance
        pdf = @(rho) exp(-log_B(k) + log(rho) * (a(k)-1) + log(1-rho) * (b(k)-1));
        
        % warn = warning('off');
        % err = (1 - quad(pdf, 0, 1, precision))
        % if err > 0.05
        %     warning('ebfret.analysis.dist.dirichlet.tau:CannotIntegrate', ...
        %         'Numerical integration of p(tau) is not accurate.');
        % end
        % E_tau(k) = quad(@(rho) -1./log(rho) .* pdf(rho), 1e-3, 1-1e-3, precision);
        % V_tau(k) = quad(@(rho) (1./log(rho)).^2 .* pdf(rho), 1e-3, 1-1e-3, precision);
        % warning(warn);
        
        % calculate most likely value
        log_pdf = @(tau) -log_B(k) - 2*log(tau) - (a(k)./tau) + log(1-exp(-1./tau)) .* (b(k)-1);
        rho = exp(linspace(log(1e-4 * E_rho(k)), log(0.5), 1e4));
        rho = [rho, 1-rho(end:-1:1)];
        tau = -1 ./ log(rho);
        [void t] = max(log_pdf(tau));
        M_tau(k) = tau(t);
    end
    if transposed
        M_tau = M_tau';
        % E_tau = E_tau';
        % V_tau = V_tau';
    end
