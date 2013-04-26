function [g_f, xi_f, bad] = jitter_filter(z_hat, g, xi, ignore)
    [T K] = size(g);
    bad = find((z_hat(3:end-1) ~= z_hat(2:end-2)) ...
               & (z_hat(1:end-3) ~= z_hat(2:end-2))) + 1;
    g_f = g;
    xi_f = xi;
    tdxs = [];
    % fix one point jumps
    if strcmp(ignore, 'spike') + strcmp(ignore, 'all')
        % time points with spikes
        tdxs = bad(z_hat(bad-1) == z_hat(bad+1));
    end
    % fix one point intermediates
    if strcmp(ignore, 'intermediate') + strcmp(ignore, 'all')
        % time points with spikes
        tdxs = cat(1, tdxs, bad(z_hat(bad-1) ~= z_hat(bad+1)));
    end
    %fprintf('bad: %d, filter: %d\n', length(bad), length(tdxs));
    % k values of preceding points
    kdxs = z_hat(tdxs-1);
    % k values of intermediate points
    ldxs = z_hat(tdxs);
    % k values of trailing points
    mdxs = z_hat(tdxs+1);
    % zero out stats for bad points
    g_f(tdxs, :) = 0; 
    xi_f(tdxs, :, :) = 0;
    % upweight everything else to ensure same number of counts
    g_f = (T-1) / (T - length(tdxs) - 1) * g_f;
    xi_f = (T-1) / (T - length(tdxs) -1) * xi_f;
    % fix transition probabilities in preceding points 
    tkk = sub2ind(size(xi_f), tdxs-1, kdxs, kdxs);
    tkl = sub2ind(size(xi_f), tdxs-1, kdxs, ldxs);
    xi_f(tkk) = xi_f(tkk) + xi_f(tkl);
    xi_f(tkl) = 0;
    % fix transition probabilities in trailing points 
    tlm = sub2ind(size(xi_f), tdxs+1, ldxs, mdxs);
    tmm = sub2ind(size(xi_f), tdxs+1, mdxs, mdxs);
    xi_f(tmm) = xi_f(tmm) + xi_f(tlm);
    xi_f(tlm) = 0;
