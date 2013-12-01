function varargout = load_sf_tracer(data_files)
    % varargout = load_sf_tracer(data_files)
    %
    % Read saved traces from SF Tracer. Each row contains a fluorphore 
    % signal for a different molecule. The first two lines contain 
    % header information, which is ignored.
    % 
    % 1:      <channel background, foreground and cutoff>
    % 2:      <column labels>
    % 3-end:  region  channel  area  length  background  intensity_1  ...  intensity_T 
    %
    %
    % Inputs
    % ------
    % 
    % data_files : [D 1] cell
    %   Filename, or filenames to load data from
    %
    % 
    % Outputs
    % -------
    %
    % channels : [N 1] cell
    %   Each output argument is a cell array containing a fluorophore 
    %   signal for two-color FRET these are [donor, acceptor]
    if isstr(data_files)
        data_files = {data_files};
    end

    for f = 1:length(data_files)
        % read sf_tracer data files 
        data = importdata(data_files{f});
        data = data.data;
        % loop over lines / fluorphores
        for l = 1:size(data)
            % molecule index is on 1st column
            % (note that numbering starts at 0)
            n = data(l, 1) + 1;
            % channel index is on 2nd column
            c = data(l, 2) + 1;
            % background intensity is on 4th column
            bg = data(l, 5);
            % signal intensities are on remaining columns
            varargout{c}{n} = data(l, 6:end)' - bg;
        end
    end
