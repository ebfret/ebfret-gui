function str = join(delim, strs)
    % str = join(delim, strs)
    %
    % Returns concatenated string using specified delimiters.
    %
    % join('/', {'some', 'directory', 'path'}) -> 'some/directory/path'
    parts = cellfun(@(s) [s, delim], strs, 'UniformOutput', false);
    str = [parts{:}];
end
