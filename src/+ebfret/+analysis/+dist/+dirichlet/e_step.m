function E_ln_theta = e_step_dir(alpha)
    % E_ln_theta = e_step_dir(alpha)
    % 
    % E-step of VBEM algorithm for Dirichlet distribution
    %
    % E[ln theta(l,k)] = Int d theta Dir(theta | alpha) ln(theta)
    %                  = psi(alpha(l,k)) - psi(sum_k alpha(l,k))
    if size(alpha,1) == prod(size(alpha))
        transposed = true;
        alpha = alpha';
    end
    E_ln_theta = bsxfun(@minus, psi(alpha), psi(sum(alpha, 2)));
    if transposed
        E_ln_theta = E_ln_theta';
    end
