function E_ln_px_z = e_step(w, x)
    % E_ln_px_z = e_step(w, x)
    % 
    % E-step of VBEM algorithm for Normal-Wishart distribution.

    % get dimensions
    [K D] = size(w.mu);

    % Expectation of log emission precision |Lambda| 
    % under q(L | w.W, w.nu) (CB 10.65, JKC 44)
    %
    % E_ln_det_L(k)  =  E[ln(|Lambda(k)|)]
    %                =  ln(|w.W|) + D ln(2) 
    %                   + Sum_d psi((w.nu(k) + 1 - d)/2)
    if D>1
        E_ln_det_L = zeros(length(w.nu), 1);  
        for k=1:length(w.nu)
          E_ln_det_L(k) = log(det(w.W(k, :, :))) + D * log(2) + ...
                          sum(psi(0.5 * (w.nu(k) + 1 - (1:D))), 2);
        end
    else
        E_ln_det_L = log(2 * w.W) + psi(0.5 * w.nu);
    end 

    % replicate w.nu w.beta w.W, and E_ln_det_L if necessary
    if length(w.beta) == 1
        w.beta = ones(K,1) * w.beta;
    end
    if length(w.nu) == 1
        w.nu = ones(K,1) * w.nu;
    end     
    if length(w.W) == 1
        w.W = bsxfun(@times, ones(K,1), w.W);
    end
    if length(E_ln_det_L) == 1
        E_ln_det_L = ones(K,1) * E_ln_det_L;
    end

    % Expectation of Mahalanobis distance Delta^2 under q(theta | w)
    % (10.64, JKC 44)
    %
    % E_Delta2(t, k) 
    %   = E[(x(t,:) - mu(k,:))' * Lambda * (x(t,:) - mu(l,:))]
    %   = D / w.beta(k) 
    %    + w.nu(k) Sum_de dx(t, d, k) W(d, e) dx(t, e, k)
    if D>1
        % dx(d, t, k) = x(t, d) - mu(k, d)
        dx = bsxfun(@minus, x', reshape(w.mu', [D 1 K]));
        % W(d, e, k) = w.W(k, d, e)
        W = permute(w.W, [2 3 1]);
        % dxWdx(t, k) = Sum_de dx(d,t,k) * W(d, e, k) * dx(e, t, k)
        dxWdx = squeeze(mtimesx(reshape(dx, [1 D T K]), ...
                                mtimesx(reshape(W, [D D 1 K]), ...
                                        reshape(dx, [D 1 T K]))));
        % note, the mtimesx function applies matrix multiplication to
        % the first two dimensions of an N-dim array, while using singleton
        % expansion to the remaining dimensions. 
        % 
        % TODO: make mtimesx usage optional? (needs compile on Linux/MacOS)
    else
        % dx(t, k) = x(t) - mu(k)
        dx = bsxfun(@minus, x, w.mu');
        % dxWdx(t, k) = Sum_de dx(t,k) * W(k) * dx(t, k)
        dxWdx = bsxfun(@times, dx, bsxfun(@times, w.W', dx));
    end
    % E_md(t, k) = D / w.beta(k) + w.nu(k) * dxWdx(t,k)
    E_Delta2 = bsxfun(@plus, (D ./ w.beta)', bsxfun(@times, w.nu', dxWdx));

    % Log expectation of p(x | z, theta) under q(theta | w)
    %
    % E_ln_px_z(t, k)
    %   = log(1 / 2 pi) * (D / 2)
    %     + 0.5 * E[ln(|Lambda(k,:,:)]]
    %     - 0.5 * E[Delta(t,k)^2]
    E_ln_px_z = log(2 * pi) * (-D / 2) ...
                + bsxfun(@minus, 0.5 * E_ln_det_L', ...
                                 0.5 * E_Delta2);
                
