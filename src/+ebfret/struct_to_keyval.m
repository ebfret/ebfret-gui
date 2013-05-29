function kv = struct_to_keyval(s, strip_empty)
    if nargin < 2
        strip_empty = false;
    end
    if ~isscalar(s) || ~isstruct(s)
        error('ebfret.struct_to_keyval:NonScalarStruct', ...
            'Argument must be a scalar struct.')
    end
    vals = struct2cell(s);
    keys = fieldnames(s);
    if strip_empty
        v = find(~[cellfun(@isempty, vals) & cellfun(@isnumeric, vals)]);
        vals = vals(v);
        keys = keys(v);
    end

    kv = cat(1, {keys{:}}, {vals{:}});
    kv = {kv{:}};
