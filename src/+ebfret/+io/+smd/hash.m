function id = hash(data, varargin)    
    import deps.ebfret.deps.datahash.datahash
    import timeseries.isvalid
    if timeseries.isvalid(data)
        id = ebfret.deps.datahash.datahash(data.data, varargin{:});
    else
        id = ebfret.deps.datahash.datahash(data, varargin{:});
    end
end
