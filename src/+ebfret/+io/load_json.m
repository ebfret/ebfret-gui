function data = load_json(filename, varargin)
    % data = load_json(filename, varargin)
    % 
    % Loads dataset written to disk in JSON format.
    %
    % Inputs
    % ------
    %
    %   filename : string
    %       Output file to read input from. If extension is .json 
    %       file assumed to be plain text.  If .json.gz input is 
    %       assumed to be compressed with gzip.
    %
    % Variable Inputs
    % ---------------
    %   
    %   gzip : false | true
    %       Gzip output. Overrides detection based on extension
    %
    %   opts : struct 
    %       Any remaining opts to be passed to the loadjson function 
    %       from the JSONlab package.
    ip = inputParser();
    ip.StructExpand = true;
    ip.KeepUnmatched = true;
    ip.addRequired('filename', @isstr);
    ip.parse(filename, varargin{:});
    args = ip.Results;
    opts = ip.Unmatched;

    % check if using gzip
    use_gzip = strcmp(filename(end-1:end), 'gz');
    if isfield(opts, 'gzip')
        use_gzip = ip.Unmatched.gzip;
        opts = rmfield(opts, 'gzip');
    end

    if use_gzip
        % read gzipped file from disk (use Java)
        import java.io.*;
        import java.util.zip.GZIPInputStream;
        buffer = BufferedReader(InputStreamReader(GZIPInputStream(FileInputStream(filename))));
        lines = {};
        l = char(buffer.readLine());
        while ~isempty(l)
            lines{end+1} = l;
            l = char(buffer.readLine());
        end
        text = strcat(lines{:});
    else
        fid = fopen(filename,'rt');
        text = fscanf(fid, '%c');
    end
    data = ebfret.deps.jsonlab.loadjson(text, opts);
