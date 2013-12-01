function smd = filter(smd, varargin)
% smd = filter(smd, varargin)
%
% Selects a subset of traces from a single-molecule dataset
%
% Variable Inputs
% ---------------
%   ids : cell array, or regular expression
%       Filters using a list of ids, or matches ids to a regular 
%       expression
%   attr : struct
%       Filters series by specified attributes. String-valued 
%       attributes are interpreted as regular expressions.
%   func : function handle
%       Custom function that returns true if trace is to be retained
%
% Jan-Willem van de Meent
% $Revision: 1.0$  $Date: 2012/10/16$

% parse inputs
ip = inputParser();
ip.StructExpand = false;
ip.addRequired('smd', @sm_isvalid);
ip.addParamValue('ids', @(i) iscell(i) | isstr(i));
ip.addParamValue('attr', struct([]), @isstruct);
ip.addParamValue('func', @(d) true, @(f) isa(f, 'function_handle'));
ip.parse(smd, varargin{:});
args = ip.Results;

data = smd.data;

% filter by id
if ~isempty(args.ids)
    if iscell(args.ids)
        [m, i] = ismember(args.ids, {data.id});
        data = [data(i)];
    elseif isstr(args.ids)
        ids = {data.id};
        msk = cellfun(@(id) any(regexp(id, args.ids)), {data.id});
        data = data(msk);
    end
end

% filter by attrs
if ~isempty(args.attr)
    if ~all(ismember(fieldnames(args.attr), fieldnames(data.attr)))
        data = data([]);
    else
        msk = ones(size(data));
        fields = fieldnames(args.attr)
        for f = 1:length(fields)
            fld = fields(f);
            if isstr(args.attr.(fld))
                % for string attributes, match regexp
                r = args.attr.(fld);
                msk = cellfun(@(a) regexp(a.(fld), r), [data.attr]);
            else
                % for other attributes do deep comparison
                msk = arrayfun(@(d) isequal(args.attr.(fld), d.attr.(fld)), data);                        
            end
            data = data(msk);
        end
    end
end

% filter by function
msk = arrayfun(@args.func, data);
data = data(msk);

% check if anything filtered
if length(data) ~= length(smd.data)
    smd.data = data;
    smd.id = ebfret.deps.datahash.datahash(data);
end


