function c = line_colors(num)
% c = line_colors(num)
%
% Returns cell array with specified number of line colors.
if num <= 6 
    c = struct2cell(ebfret.plot.named_colors());
    c = {c{1:num}};
else
    c = mat2cell(ebfret.plot.spectrum(num), ones(num,1), 3);
end
c = {c{:}};