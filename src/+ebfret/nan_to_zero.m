function A = nan2zero(A)
    A(isnan(A)) = 0;