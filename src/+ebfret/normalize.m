function [A0, C] = normalize(A, dim)
% [A0, C] = normalize(A, dim)
%
% Returns normalized array A0 and normalization const C. Array is
% along single axis if dim is specified, or along all axes if not.
if nargin < 2
  C = sum(A(:));
  A0 = A ./ (C + (C==0));  
else
  C = sum(A, dim);
  A0 = bsxfun(@rdivide, A, C + (C==0));
end
