function to_csv(report, filename, sep, tab)
    if nargin < 2
        [fname, fpath, findex] = ...
            uiputfile({'*.csv', 'ebFRET analysis summary (.csv)';});
        filename = sprintf('%s/%s', fpath, fname);
    end
    if nargin < 3
        sep = ',';
    end
    if nargin < 4
        tab = '    ';
    end
    function field = indent(field, level)
        for l = 1:level
            field = sprintf('%s%s', tab, field);
        end
    end
    function lines = parse(node, level)
        if nargin < 2
            level = 0;
        end
        fields = fieldnames(node);
        lines = {};
        for f = 1:length(fields)
            r0 = length(lines);
            lines{r0+1,1} = {indent(fields{f}, level)};
            if isstruct(node(1).(fields{f}))
                lines = cat(1, lines, parse([node.(fields{f})], level+1));
            elseif isstr(node(1).(fields{f}))
                lines{r0+1,1}{2} = sprintf([sep, '%s'], node.(fields{f}));
            elseif iscell(node(1).(fields{f}))
                for n = 1:length(node)
                    lines{r0+1,1}{1+n} = sprintf([sep, '%s'], node(n).(fields{f}){:});
                end
            elseif isvector(node(1).(fields{f}))
                for n = 1:length(node)
                    lines{r0+1,1}{1+n} = sprintf([sep, '%.3e'], node(n).(fields{f}));
                end
            elseif ismatrix(node(1).(fields{f}))
                for n = 1:length(node)
                    data = node(n).(fields{f});
                    for r = 1:size(data, 1);
                        lines{r0+r,1}{1+n} = sprintf([sep, '%.3e'], data(r,:));    
                    end
                end
            end
        end 
    end 
    lines = parse(report);
    fid = fopen(filename, 'wt');
    for l = 1:length(lines)
        fprintf(fid, sprintf('%s\n', cat(2, lines{l}{:})));
    end
    fclose(fid);
end