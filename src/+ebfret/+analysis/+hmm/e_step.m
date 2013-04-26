function [E_ln_pi, E_ln_A, E_ln_px_z] = e_step(w, x)
    % [E_ln_pi, E_ln_A, E_p_x_z] = e_step(w)
    % 
    % E-step of VBEM algorithm for a HMM with gaussian-distributed
    % observations.

    % get dimensions
    [K D] = size(w.mu);

    % Expectation of log intial state priors pi under q(pi | w.pi) 
    % (MJB 3.69, CB 10.66, JCK 41)
    %
    % E[ln(w.pi(k))]  =  Int d pi  Dir(pi | w.pi) ln(pi)
    %                 =  psi(w.pi(k)) - psi(Sum_l w.pi(l)))
    E_ln_pi = psi(w.pi) - psi(sum(w.pi)); 

    % Expectation of log transition matrix A under q(A | w.A) 
    % (MJB 3.70, JCK 42)
    %
    % E_ln_A(k, l)  =  psi(w.A(k,l)) - psi(Sum_l w.A(k,l))
    E_ln_A = bsxfun(@minus, psi(w.A), psi(sum(w.A, 2)));

    % Expectation of log precision and log emission probability
    E_ln_px_z = ebfret.analysis.dist.normwish.e_step(w, x);
