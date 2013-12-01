function write_report(filename, report, varargin)
    ip = inputParser();
    ip.StructExpand = true;
    ip.addParamValue('separator', ',', @isstr);       
    ip.addParamValue('indent', '    ', @isstr);       
    ip.parse(varargin{:});
    args = ip.Results;

    if nargin < 3
        args.separator = ',';
    end
    if nargin < 4
        args.indent = '    ';
    end
    function field = indent(field, level)
        for l = 1:level
            field = sprintf('%s%s', args.indent, field);
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
                lines{r0+1,1}{2} = sprintf([args.separator, '%s'], node.(fields{f}));
            elseif iscell(node(1).(fields{f}))
                for n = 1:length(node)
                    lines{r0+1,1}{1+n} = sprintf([args.separator, '%s'], node(n).(fields{f}){:});
                end
            elseif isvector(node(1).(fields{f}))
                for n = 1:length(node)
                    lines{r0+1,1}{1+n} = sprintf([args.separator, '%.3e'], node(n).(fields{f}));
                end
            elseif ismatrix(node(1).(fields{f}))
                for n = 1:length(node)
                    data = node(n).(fields{f});
                    for r = 1:size(data, 1);
                        lines{r0+r,1}{1+n} = sprintf([args.separator, '%.3e'], data(r,:));    
                    end
                end
            end
        end 
    end 
    function num = count_sep(elem)
        if isstr(elem)
            num = length(strfind(elem, args.separator));
        else
            num = 0;
        end
    end
    function lines = align_fields(lines)
        for l = 2:length(lines)
            for c = 2:length(lines{l})
                if ~count_sep(lines{l}{c})
                    n = count_sep(lines{l-1}{c});
                    lines{l}{c} = repmat(args.separator, [1 n]);
                end
            end
        end
    end
    lines = align_fields(parse(report));
    fid = fopen(filename, 'wt');
    for l = 1:length(lines)
        fprintf(fid, sprintf('%s\n', cat(2, lines{l}{:})));
    end
    fclose(fid);
end