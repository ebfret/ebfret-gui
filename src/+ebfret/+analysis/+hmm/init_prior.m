function u = init_prior(theta, counts)
	% u = init_prior(theta, counts)
	%
    % Returns prior parameters for an HMM with normally 
    % distributed observations based on a set of expected
    % values and prior counts.
    %
    %
    % Inputs
    % ------
    % 
    % theta : struct
    %   Expected parameter values under prior
    %   .mu : [K 1]
    %       Means of observation distribution
    %   .lambda : [K 1]
    %       Precisions of observation distribution
    %   .tau : [K 1]
    %       State dwell time
    %
    % counts : struct
    %   Prior strength (in number of of observations)
    %   for each component.
    %
    % Outputs
    % -------
    %
    % u : struct
    %   Prior parameters 

    % construct expectation value for transition matrix from 
    % expected dwell times
    K = length(theta.tau);
    if K > 1
        Akk = exp(-1 ./ theta.tau);
        E_A = bsxfun(@times, (1-eye(K)), (1-Akk(:)) / (K-1))  ...
              + bsxfun(@times, eye(K), Akk(:));

        % get expected initial probabilities from eigenvalue of 
        % transition matrix
        [v, l] = eig(E_A');
        [void k] = max(diag(l));
        E_pi = ebfret.normalize(v(:,k));
    else
        E_A = 1;
        E_pi = 1;
    end
        
    % construct prior
    u.pi = mean(counts.tau(:)) .* E_pi(:);
    u.A = bsxfun(@times, counts.tau(:), E_A);
    u.mu = theta.mu(:);
    u.beta = counts.mu(:);
    u.W = theta.lambda(:) ./ counts.lambda(:);
    u.nu = counts.lambda(:);
