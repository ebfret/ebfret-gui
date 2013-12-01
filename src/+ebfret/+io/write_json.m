function write_json(filename, data, varargin)
    % write_json(filename, data, varargin)
    % 
    % Writes dataset to disk in JSON format.
    %
    % Inputs
    % ------
    %
    %   filename : string
    %       Output file to write to. If extension is .json file is
    %       written to plain text. If .json.gz output is compressed
    %       with gzip before writing to disk.
    %
    %   data : struct
    %       Data to write to disk.
    %
    % Variable Inputs
    % ---------------
    %   
    %   gzip : false | true
    %       Gzip output. Overrides detection based on extension
    %
    %   opts : struct 
    %       Any remaining opts to be passed to the savejson function 
    %       from the JSONlab package.
    ip = inputParser();
    ip.StructExpand = true;
    ip.KeepUnmatched = true;
    ip.addRequired('filename', @isstr);
    ip.addRequired('data', @isstruct);
    ip.parse(filename, data, varargin{:});
    args = ip.Results;
    opts = ip.Unmatched;

    % check if using gzip
    use_gzip = strcmp(filename(end-1:end), 'gz');
    if isfield(opts, 'gzip')
        use_gzip = ip.Unmatched.gzip;
        opts = rmfield(opts, 'gzip');
    end

    % convert data to json
    json = ebfret.deps.jsonlab.savejson('', args.data, opts);
    
    % write to disk (use Java)
    import java.io.*;
    out = FileOutputStream(args.filename);
    if use_gzip
        import java.util.zip.GZIPOutputStream;
        writer = OutputStreamWriter(GZIPOutputStream(out));
    else
        writer = OutputStreamWriter(out);
    end
    writer.write(json);
    writer.close();
    out.close();