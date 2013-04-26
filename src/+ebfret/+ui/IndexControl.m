classdef IndexControl < hgsetget
    properties
        min;
        max;
        value;
        handles;
        parent_callback;
    end
    methods
        function self = IndexControl(varargin)
            % self = indexControl(varargin)
            %
            % Variable Arguments
            % ------------------
            %   'Parent' : handle
            %   'Position' : [l b h w]
            %   'Units' : {'pixels', 'normalized'}
            %   'Callback' : {@function, arguments}
            %   'Min' : int
            %   'Max' : int
            %   'Value' : int

            ip = inputParser();
            ip.StructExpand = true;
            ip.addParamValue('parent', gco(), @isnumeric);       
            ip.addParamValue('position', [0 0 1 1], @isnumeric);       
            ip.addParamValue('units', 'Normalized', ...
                             @(s) any(strcmpi(s, ...
                                {'inches', 'centimeters', 'points', ...
                                 'normalized', 'pixels'})));       
            ip.addParamValue('parent_callback', @(n) show_series(gco(), n));
            ip.addParamValue('min', 0, @isscalar);
            ip.addParamValue('max', 0, @isscalar);
            ip.addParamValue('value', 0, @isscalar);
            ip.parse(varargin{:});
            args = ip.Results;

            % round limits to nearest integer
            args.min = round(args.min);
            args.max = round(args.max);

            % initialize enclosing frame
            handles.panel ...
                = uipanel('parent', args.parent, ...
                          'backGroundColor', [0.95 0.95 0.95], ...
                          'position', args.position, ...
                          'units', args.units);

            % horizontal and vertical padding, slider height (normalized units)
            pos = getpixelposition(handles.panel);
            hp = 4 / pos(3);
            vp = 4 / pos(4);
            sh = 22 / pos(4);

            % intialize control elements in frame: 
            % text edit, slider, prev & next buttons
            handles.indexEdit ...
                = uicontrol('parent', handles.panel, ... 
                            'style', 'edit', ...
                            'backGroundColor', [1 1 1], ...
                            'callback', @(source, event) callback(self, source, event), ...
                            'units', 'normalized', ...
                            'position', [0.2*hp, vp, 0.1-0.2*hp, 1-2*vp]);
            handles.slider ...
                = uicontrol('parent', handles.panel, ... 
                            'style', 'slider', ...
                            'min', args.min-eps, ...
                            'max', args.max+eps, ...
                            'value', args.value, ...
                            'callback', @(source, event) callback(self, source, event), ...
                            'units', 'normalized', ...
                            'position', [0.1+hp, vp, 0.7-hp, sh]);
            handles.prevButton ...
                = uicontrol('parent', handles.panel, ...
                            'style', 'push', ...
                            'string', 'Prev', ...
                            'callback', @(source, event) callback(self, source, event), ...
                            'units', 'normalized', ...
                            'position', [0.8+hp, vp, 0.1-hp, 1-2*vp]);
            handles.nextButton ...
                = uicontrol('parent', handles.panel, ...
                            'style', 'push', ...
                            'string', 'Next', ...
                            'callback', @(source, event) callback(self, source, event), ...
                            'units', 'normalized', ...
                            'position', [0.9+hp, vp, 0.1-2*hp, 1-2*vp]);
            set(self, 'handles', handles);
            set(self, 'parent_callback', args.parent_callback);

            set(self, 'max', nan, 'min', nan, 'value', nan);
            set_limits(self, args.max, args.min);
            set_value(self, args.value);
        end
        % callback wrapper
        function callback(self, source, event)
            handles = get(self, 'handles');
            state.value = get(self, 'value');
            switch source
                % get new value
                case handles.indexEdit
                    value = round(str2num(get(handles.indexEdit, 'string')));
                case handles.slider
                    value = round(get(handles.slider, 'value'));
                case handles.prevButton
                    value = state.value - 1;
                case handles.nextButton
                    value = state.value + 1;
            end
            % update widget value
            set_value(self, value);
        end
        function set_value(self, value)
            % ensure value in range
            value = max(min(value, self.max), self.min);
            % check if update needed
            if value ~= self.value;
                % update edit box and slider value
                set(self.handles.indexEdit, 'string', sprintf('%d', value));
                set(self.handles.slider, 'value', value);
                % disable / enable next button if at max / previously at max
                if (value == self.max) | isnan(value)
                    set(self.handles.nextButton, 'enable', 'off');
                end
                if (self.value == self.max) & (value < self.max)
                    set(self.handles.nextButton, 'enable', 'on');
                end
                % disable / enable prev button if at min / previously at min
                if (value == self.min) | isnan(value)
                    set(self.handles.prevButton, 'enable', 'off');
                end
                if (self.value == self.min) & (value < self.min)
                    set(self.handles.prevButton, 'enable', 'on');
                end
                % set new state value
                set(self, 'value', value);
                % send update signal to main window
                parent_callback = get(self, 'parent_callback');
                parent_callback(value);
            end
        end
        function set_limits(self, max_val, min_val)
            if nargin < 3
                min_val = 1;
            end
            max_val = round(max_val);
            min_val = round(min_val);
            % update slider limits
            set(self.handles.slider, ...
                'min', min_val-eps, 'max', max_val+eps);
            % enable / disable controls as needed
            if max_val > min_val
                set(self.handles.indexEdit, 'enable', 'on'); 
                set(self.handles.slider, 'enable', 'on'); 
                set(self.handles.prevButton, 'enable', 'on'); 
                set(self.handles.nextButton, 'enable', 'on'); 
            else
                set(self.handles.indexEdit, 'enable', 'off'); 
                set(self.handles.slider, 'enable', 'off'); 
                set(self.handles.prevButton, 'enable', 'off'); 
                set(self.handles.nextButton, 'enable', 'off'); 
            end
            % update value
            set_value(self, self.value);
        end
    end
end