function smd = merge(varargin)
    % smd = filter(varargin)
    %
    % Merges supplied single-molecule datasets into a single 
    % structure. Any dataset-level attributes which are not 
    % identical for all datasets are made attributes of individual 
    % time series. Conversely, all series level attributes that are 
    % identical for the entire dataset, are made dataset level 
    % attributes.
    %
    % Jan-Willem van de Meent
    % $Revision: 1.0$  $Date: 2012/10/16$

    smds = cellfun(@(smd) smd(:), varargin, 'UniformOutput', false);
    smds = cat(1, smds{:});

    % check that data types are the same
    if ~all(strcmp(smds(1).type, {smds.type}))
        error('SMD:TypeMismatch', ...
              'Data types are not consistent for all arguments.')
    end    

    % check that label names are the same
    for d = 2:length(smds)
        if ~all(strcmp(smds(1).columns, smds(d).columns))
            error('SMD:ColumnsMismatch', ...
                  'Column labels are not consistent for all arguments.')
        end
    end

    % copy dataset level attrs to individual traces
    for d = 1:length(smds)
        attrs = [smds(d).data.attr];
        for f = fieldnames(smds(d).attr)'
            [attrs.(char(f))] = deal(smds(d).attr.(char(f)));
        end
        for n = 1:length(smds(d).data)
            smds(d).data(n).attr = attrs(n);
        end
    end

    % initialize merged dataset
    smd = struct('id', {}, 'attr', {}, 'columns', {}, 'data', {});
    smd(1).columns = smds(1).columns;
    smd(1).data = struct('id', {}, 'attr', {}, 'index', {}, 'values', {});

    mmin = 0;
    for d = 1:length(smds)
        for n = 1:length(smds(d).data)
            m = mmin + n;
            [smd.data(m).id] = smds(d).data(n).id;        
            [smd.data(m).index] = smds(d).data(n).index;        
            [smd.data(m).values] = smds(d).data(n).values;
            for f = fieldnames(smds(d).data(1).attr)'
                [smd.data(m).attr.(char(f))] = ...
                    smds(d).data(n).attr.(char(f));
            end
        end
        mmin = mmin + length(smds(d).data);
    end

    % move all identical attributes back to dataset level
    attrs = [smd.data.attr];
    for f = fieldnames(attrs)'
        f = char(f);
        if all(arrayfun(@(a) isequal(a.(f), attrs(1).(f)), attrs))
            smd.attr.(f) = attrs(1).(f);
            attrs = rmfield(attrs, f);
        end
    end
    for n = 1:length(smd.data)
        smd.data(n).attr = attrs(n);
    end
    
    % calculate id from hash
    smd.id = ebfret.deps.datahash.datahash(smd.data);
end