function lines = cat(dim, varargin)
    % lines = cat(dim, a, b)
    % 
    % Concatenates sets of lines along specified dimension
    lines = struct([]);
    for a = 1:length(varargin)
        arg = varargin{a};
        fn = fieldnames(arg);
        rng = (length(lines)+1):(length(lines)+length(arg));
        for f = 1:length(fn)
            [lines(rng).(fn{f})] = arg.(fn{f});
        end
    end
