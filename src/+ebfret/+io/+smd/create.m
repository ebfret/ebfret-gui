function smd = create(data, columns, varargin)
% smd = create(data, columns, varargin)
%
% Creates a single-molecule dataset structure from supplied data.
%
% Inputs
% ------
%   data : N x 1 cell,  M x D numeric,  N x T x D numeric
%       Time series data. May be formatted in 3 ways
%       1. As a cell array, where each data{n} is a T{n} x D array
%       2. As a T x D flat array, e.g. of form [id, time, values]
%       3. As a N x T x D array 
%   columns : D x 1 cell   
%       Column labels for time series data, eg {'id', 'time', 'fret'}. 
%       The following labels have special meanings:
%       'index' : used to specify smd.data(n).index
%       'id' : used to specify smd.data(n).id
%
% Variable Inputs
% ---------------
%   index : N x T  or N x 1 cell of T{n} x 1 float
%       May be used to separately specify time index values. Overrides 
%       entries supplied in data argument.
%   type : string
%       Type id for data. Constructed from columns if left blank.
%   id : string
%       Unique id for dataset. Computed using hash if not specified.
%   attr : struct
%       May be used to provide additional annotation of data.
%   data_ids : N x 1 cell 
%       Unique ids for traces. Computed using hash if not specified.
%   data_attrs : N x 1 struct 
%       May be used to store additional information specific to each
%       individual time series.
%
% Jan-Willem van de Meent
% $Revision: 1.0$  $Date: 2012/10/16$

% parse inputs
ip = inputParser();
ip.StructExpand = false;
ip.addRequired('data');
ip.addRequired('columns', @iscell);
ip.addParamValue('index', {});
ip.addParamValue('type', '', @isstr);
ip.addParamValue('id', '', @isstr);
ip.addParamValue('attr', struct(), @isstruct);
ip.addParamValue('data_ids', {}, @iscell);
ip.addParamValue('data_attrs', struct(), @isstruct);
ip.parse(data, columns, varargin{:});
args = ip.Results;

% parse columns argument
[m, i] = ismember('index', lower(args.columns));
index_label = m * i;
[m, i] = ismember('id', lower(args.columns));
id_label = m * i;
value_columns = setdiff(1:length(args.columns), [index_label, id_label]);
D = length(value_columns);

% helper function for parsing data argument
function d = parse_data(data, dfs, tf, idf)
    if idf
        % split by id if necessary
        id = data(:, idf);
        ids = unique(id);
        for n = 1:length(ids)
            d(n).id = num2str(ids(n));
            i = find(id == ids(n));
            if tf 
                d(n).index = data(i, tf);
            else
                d(n).index = (1:length(i))';
            end
            d(n).values = data(i, dfs);
        end
    else
        % assume single trace with blank id
        d.id = '';        
        if tf
            d.index = data(:, tf);
        else
            d.index = (1:length(data))';
        end
        d.values = data(:, dfs);
    end
end

% parse data argument
if iscell(args.data)
    data = cellfun(@(d) parse_data(d, value_columns, ...
                                   index_label, id_label), ... 
                   {args.data{:}});
elseif isnumeric(args.data)
    switch ndims(args.data)
        case 2
            data = parse_data(args.data, value_columns, index_label, id_label);
        case 3
            data = ...
                arrayfun(@(n) parse_data(squeeze(args.data(n,:,:)), ...
                                         value_columns, index_label, id_label), ...
                         1:size(args.data, 1));
        otherwise
            error('SMD:InvalidInput', ...
                  'The data argmument must have either have size [N T D] or [M D]');
    end
else
    error('SMD:InvalidInput', ...
          'The "data" argument must either be a cell or numeric array.');
end

% set indexes if specified
if not(isempty(args.index))
    if iscell(args.index)
        index = args.index;
    elseif isnumeric(args.index)
        index = num2cell(args.index, 2);
    else
        error('SMD:InputsNotAligned', ...
              'The "index" argument must either be a cell or numeric array.');
    end
    if length(args.index) ~= length(data)
        error('SMD:InvalidInput', ...
              'Number of specified index values (%d) does not match number of series parsed from data (%d)', ...
              length(index), length(data));
    end
    for n = 1:length(data)
        data(n).index = index{n}(:);
    end
end

% set id's if specified
if not(isempty(args.data_ids))
    if length(args.data_ids) ~= length(data)
        error('SMD:InputsNotAligned', ...
              'Number of specified data_ids (%d) does not match number of time series parsed from data (%d)', ...
              length(args.data_ids), length(data));
    end
    for n = 1:length(data)
        data(n).id = args.data_ids{n};
    end
end

% set attr's if specified
if length(args.data_attrs) > 1
    if length(args.data_attrs) ~= length(data)
        error('SMD:InputsNotAligned', ...
              'Number of specified data_attrs (%d) does not match number of series parsed from data (%d)', ...
              length(args.data_attrs), length(data));
    end
    for n = 1:length(data)
        data(n).attr = args.data_attrs(n);
    end
else 
    [data.attr] = deal(struct());
end

% replace any empty time series ids with hashes
for n = 1:length(data)
    if isempty(data(n).id)
        data(n).id = ebfret.deps.datahash.datahash(data(n).values);
    end
end

% calculate set id from hash if necessary
if isempty(args.id)
    id = ebfret.deps.datahash.datahash(data);
else
    id = args.id;
end

% type contains column labels if unspecified
if isempty(args.type)
    args.type = [sprintf('%s-', columns{1:end-1}), columns{end}];
end

% assign data structure
smd = struct();
smd.type = args.type;
smd.id = id;
smd.attr = args.attr;
smd.columns = args.columns(value_columns);
smd.data = data;
end