function [donor acceptor labels] = load(data_files, varargin)
% Loads traces from a number of data files into a single dataset.
%
% Inputs
% ------
%
% data_files (1xD cell)
%   Names of datasets to analyse (e.g. {'file1.dat' 'file2.dat'}).
%   Datasets may have one of two possible formats.
%
%   dat : (Tx3 numeric)
%       dat is formatted as a stacked table with series index on 
%       the first, donor signal on the second, and acceptor
%       signal on the 3rd column
%
%   dat : (Tx2N numeric)
%       dat is formatted as a matrix with donor signals on odd columns
%       and acceptor signals on even colums
%
%   dat : (Nx1 cell)
%       Each trace dat{n} is a T(n)x2 matrix with donor on the first
%       and acceptor on the second column.
%
% Variable Inputs
% ---------------
%
% 'variable' (string, default:empty)
%   If specified, assume data is saved as a matlab file,
%   with data stored under the specified variable name.
%
% 'has_labels' (boolean, default:true)
%   Assume first row contains trace labels
%
% 'strip_first' (boolean, default:false)
%   Discard first time point in each trace
%
%
% Outputs
% -------
%
% donor : [1 N] cell
%   Donor fluorophore signals
%
% acceptor : [1 N] cell
%   Acceptor fluorophore signals
%
% labels : [1 N] cell
%   Labels for each time series
%
% Jan-Willem van de Meent
% $Revision: 1.10 $  $Date: 2012/02/10$

% parse inputs
ip = inputParser();
ip.StructExpand = true;
ip.addRequired('data_files', @(d) iscell(d) | isstr(d));
ip.addParamValue('variable', '', @isstr);
ip.addParamValue('has_labels', true, @isscalar);
ip.addParamValue('strip_first', 0, @isscalar);
ip.parse(data_files, varargin{:});

% collect inputs
args = ip.Results;
if isstr(args.data_files)
    data_files = {args.data_files};
else
    data_files = args.data_files;
end

for d = 1:length(data_files)
    % load dataset
    if isempty(args.variable)
        dat = load(data_files{d});
    else
        dat = load(data_files{d}, args.variable);
        dat = dat.(args.variable);
    end

    % check if data is in matrix format
    if isnumeric(dat)
        % convert data into cell array
        if size(dat, 2) == 3
            ns = dat(:,1);
            for n = unique(ns(:)')
                raw_data{n} = dat(find(n == ns), 2:3);
            end
        else 
            raw_data = mat2cell(dat, size(dat,1), 2 * ones(size(dat,2) / 2, 1));
        end
    else
        % assume a cell array with dat{n} = [donor, acceptor]
        raw_data = dat;
    end

    for n = 1:length(raw_data)
        % get trace
        rw = raw_data{n};

        if ~isempty(rw)
            % strip labels if necessary
            if args.has_labels
                labels(n,1) = rw(1, 1);
                rw = rw(2:end, :);
            else
                labels(n,1) = n;
            end

            % strip first point (often bad data)
            if args.strip_first
                rw = rw(2:end, :);
            end

            % store trace again
            raw_data{n} = rw;
        else
            labels(n,1) = nan;
        end
    end

    % strip empty traces
    n = find(~cellfun(@isempty, raw_data));
    raw_data = raw_data(n);
    labels = labels(n);

    % store raw data in output
    data(d).donor = cellfun(@(rd) rd(:,1), raw_data, 'UniformOutput', false);
    data(d).acceptor = cellfun(@(rd) rd(:,2), raw_data, 'UniformOutput', false);
    data(d).labels = labels;
end

% concatenate data
donor = cat(1, data.donor);
acceptor = cat(1, data.acceptor);
labels = cat(1, data.labels);