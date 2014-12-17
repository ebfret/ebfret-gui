classdef PlotPanel < hgsetget
    properties
        handles
        plots
    end
    methods
        function self = PlotPanel(varargin)
            % self = PlotPanel(varargin)
            %
            % Variable Arguments
            % ------------------
            %   'parent' : handle
            %   'position' : [l b h w]
            %   'units' : {'pixels', 'normalized'}
            %   'title' : string
            %   'axes' : {'name', [l b h w], ...}
            ip = inputParser();
            ip.StructExpand = true;
            ip.addParamValue('parent', gco(), @ishandle);       
            ip.addParamValue('position', [0 0 1 1], @isnumeric);       
            ip.addParamValue('title', 'Time Series', @isstr);       
            ip.addParamValue('units', 'Normalized', ...
                             @(s) any(strcmpi(s, ...
                                {'inches', 'centimeters', 'points', ...
                                 'normalized', 'pixels'})));       
            ip.addParamValue('axes', {});       
            ip.addParamValue('properties', {});       
            ip.parse(varargin{:});
            args = ip.Results;

            % create panel
            handles.panel ...
                = uipanel(args.parent, ...
                          'title', args.title, ...
                          'position', args.position);

            % create axes
            if iscell(args.properties)
                properties = struct(args.properties{:});
            else
                properties = args.properties;
            end
            if iscell(args.axes)
                args.axes = struct(args.axes{:});
            end
            fields = fieldnames(args.axes);
            if isscalar(properties)
                properties(1:length(fields)) = properties;
            end
            for f = 1:length(fields)
                fld = fields{f};
                self.handles.axes.(fld) = ...
                    axes('parent', handles.panel, ...
                         'outerposition', args.axes.(fld));
                % set axes properties
                lprops = {'xlabel', 'ylabel', 'zlabel', 'title'};
                props = fieldnames(properties);
                for p = 1:length(props)
                    if any(strcmpi(props{p}, lprops))
                        t = get(self.handles.axes.(fld), props{p});
                        set(t, 'string', properties(f).(props{p}));
                        fprops = intersect(props, {'fontangle', 'fontname', 'fontunits', 'fontsize', 'fontweight'});
                        for q = 1:length(fprops)
                            set(t, fprops{q}, properties(f).(fprops{q}));
                        end 
                    else
                        set(self.handles.axes.(fld), props{p}, properties(f).(props{p}));
                    end
                end            
            end
            
            % initialize (empty) plot data
            clear_plots(self, fieldnames(self.handles.axes));
        end
        % updates plot data
        function set_plots(self, varargin)
            args = struct(varargin{:});
            fields = fieldnames(args);
            plots = get(self, 'plots');
            for f = 1:length(fields)
                if iscell(args.(fields{f}))
                    args.(fields{f}) = struct(args.(fields{f}){:});
                end
                plots.(fields{f}) = args.(fields{f});
            end
            set(self, 'plots', plots);
            refresh(self, fields);
        end
        % updates axis properties (see axes)
        function set_props(self, varargin)
            args = struct(varargin{:});
            axes = fieldnames(args);
            handles = get(self, 'handles');
            for a = 1:length(axes)
                if iscell(args.(axes{a}))
                    args.(axes{a}) = struct(args.(axes{a}){:});
                end
                set(handles.axes.(axes{a}), args.(axes{a}));
            end
            % set(self, 'props', props);
            % refresh(self, axes);
        end
        % clears plot data
        function clear_plots(self, fields)
            get(self, 'plots');
            if nargin < 2 
                fields = fieldnames(self.handles.axes);
            elseif isstr(fields)
                fields = {fields};
            end
            for f = 1:length(fields)
                plots.(fields{f}) = struct([]);
            end
            set(self, 'plots', plots);
            refresh(self, fields);
        end
        % replots panels
        function refresh(self, fields)
            plots = get(self, 'plots');
            handles = get(self, 'handles');
            if nargin < 2
                fields = fieldnames(handles.axes);
            elseif isstr(fields)
                fields = {fields};
            end
            for f = 1:length(fields)
                % clear axis
                s = warning('off');
                cla(handles.axes.(fields{f}));
                warning(s);
                % do plots
                for p = plots.(fields{f})(:)';
                    p.parent = handles.axes.(fields{f});
                    keyvals = ebfret.struct_to_keyval(p, true);
                    line(keyvals{:});
                end
                % % set axis properties
                % props = getfield(get(self, 'props'), fields{f});
                % props = cat(1, fieldnames(props)', struct2cell(props)');
                % if ~isempty(props)
                %     set(handles.axes.(fields{f}), props{:});
                % end
            end
        end
    end
end
