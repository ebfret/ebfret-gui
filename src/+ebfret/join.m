function str = join(delim, strs, last)
    % str = join(delim, strs, last)
    %
    % Returns concatenated string using specified delimiters.
    %
    % join('/', {'some', 'directory', 'path'}) -> 'some/directory/path'
    if nargin < 3
    	last = false;
    end
    parts = cellfun(@(s) {s, delim}, strs, 'UniformOutput', false);
    parts = cat(2, parts{:});
    str = [parts{1:end-1}];
end
