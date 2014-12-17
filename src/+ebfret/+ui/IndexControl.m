classdef IndexControl < hgsetget
    properties
        min;
        max;
        value;
        handles;
        callback;
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
            ip.addParamValue('parent', gco(), @ishandle);       
            ip.addParamValue('position', [0 0 1 1], @isnumeric);       
            ip.addParamValue('units', 'Normalized', ...
                             @(s) any(strcmpi(s, ...
                                {'inches', 'centimeters', 'points', ...
                                 'normalized', 'pixels'})));       
            ip.addParamValue('callback', @(value) nan);
            ip.addParamValue('min', 0, @isscalar);
            ip.addParamValue('max', 0, @isscalar);
            ip.addParamValue('value', 0, @isscalar);
            ip.addParamValue('title', '', @isstr);
            ip.parse(varargin{:});
            args = ip.Results;

            % round limits and value to nearest integer
            args.min = round(args.min);
            args.max = round(args.max);
            args.value = round(args.value);

            % initialize enclosing frame
            handles.panel ...
                = uipanel('parent', args.parent, ...
                          'backGroundColor', [0.95 0.95 0.95], ...
                          'position', args.position, ...
                          'units', args.units, ...
                          'title', args.title);

            % horizontal and vertical padding, slider height (normalized units)
            pos = getpixelposition(handles.panel);
            hp = 4 / pos(3);
            vp = 1 / pos(4);
            sh = 28 / pos(4);

            % intialize control elements in frame: 
            % text edit, slider, prev & next buttons
            self.handles.indexEdit ...
                = uicontrol('parent', handles.panel, ... 
                            'style', 'edit', ...
                            'backGroundColor', [1 1 1], ...
                            'callback', @(source, event) self.control_callback(source, event), ...
                            'units', 'normalized', ...
                            'position', [0.2*hp, vp, 0.1-0.2*hp, 1-2*vp]);
            self.handles.slider ...
                = uicontrol('parent', handles.panel, ... 
                            'style', 'slider', ...
                            'min', args.min-eps, ...
                            'max', args.max+eps, ...
                            'value', args.value, ...
                            'sliderstep', [0 0.1], ...
                            'callback', @(source, event) self.control_callback(source, event), ...
                            'units', 'normalized', ...
                            'position', [0.1+hp, vp, 0.9-2*hp, sh]);
            % self.handles.prevButton ...
            %     = uicontrol('parent', handles.panel, ...
            %                 'style', 'push', ...
            %                 'string', 'Prev', ...
            %                 'callback', @(source, event) self.control_callback(source, event), ...
            %                 'units', 'normalized', ...
            %                 'position', [0.8+hp, vp, 0.1-hp, 1-2*vp]);
            % self.handles.nextButton ...
            %     = uicontrol('parent', handles.panel, ...
            %                 'style', 'push', ...
            %                 'string', 'Next', ...
            %                 'callback', @(source, event) self.control_callback(source, event), ...
            %                 'units', 'normalized', ...
            %                 'position', [0.9+hp, vp, 0.1-2*hp, 1-2*vp]);
            self.set_prop('min', args.min, 'max', args.max, 'value', args.value, 'callback', args.callback);
        end
        % callback wrapper
        function control_callback(self, source, event)
            switch source
                % get new value
                case self.handles.indexEdit
                    value = round(str2num(get(self.handles.indexEdit, 'string')));
                case self.handles.slider
                    value = round(get(self.handles.slider, 'value'));
                % case self.handles.prevButton
                %     value = self.value - 1;
                % case self.handles.nextButton
                %     value = self.value + 1;
            end
            if value ~= self.value
                % update widget value
                self.set_prop('value', value);
                % run callback
                self.callback(value);
            end
        end
        % function value = set_value(self, value)
        %     % ensure value in range
        %     value = max(min(value, self.max), self.min);
        %     % update edit box and slider value
        %     set(self.handles.indexEdit, 'string', sprintf('%d', value));
        %     set(self.handles.slider, 'value', value);
        %     % check if update needed
        %     if value ~= self.value;
        %         % disable / enable next button if at max / previously at max
        %         if (value == self.max) | isnan(value)
        %             set(self.handles.nextButton, 'enable', 'off');
        %         end
        %         if (value < self.max)
        %             set(self.handles.nextButton, 'enable', 'on');
        %         end
        %         % disable / enable prev button if at min / previously at min
        %         if (value == self.min) | isnan(value)
        %             set(self.handles.prevButton, 'enable', 'off');
        %         end
        %         if (value > self.min)
        %             set(self.handles.prevButton, 'enable', 'on');
        %         end
        %         % set new state value
        %         set(self, 'value', value);
        %         % send update signal to main window
        %         parent_callback = get(self, 'parent_callback');
        %         parent_callback(value);
        %     end
        % end
        function [value, lim] = set_prop(self, varargin)
            properties = struct(varargin{:});
            if isfield(properties, 'callback')
                self.callback = properties.callback;
            end
            if isfield(properties, 'min')
                self.min  = round(properties.min);
                set(self.handles.slider, 'min', self.min-eps);
                if self.max > self.min
                    set(self.handles.slider, 'sliderstep', [1./(self.max-self.min) 1./(self.max-self.min)])
                end
                self.set_prop('value', max(self.value, self.min));
            end
            if isfield(properties, 'max')
                self.max  = round(properties.max);
                set(self.handles.slider, 'max', self.max+eps);
                if self.max > self.min
                    set(self.handles.slider, 'sliderstep', [1./(self.max-self.min) 1./(self.max-self.min)])
                end
                self.set_prop('value', min(self.value, self.max));
            end
            if isfield(properties, 'value')
                self.value = max(min(round(properties.value), self.max), self.min);
                set(self.handles.indexEdit, 'string', sprintf('%d', self.value));
                set(self.handles.slider, 'value', self.value);
            end
            % enable / disable controls as needed
            if self.max == self.min
                set(self.handles.indexEdit, 'enable', 'off'); 
                set(self.handles.slider, 'enable', 'off'); 
                % set(self.handles.prevButton, 'enable', 'off'); 
                % set(self.handles.nextButton, 'enable', 'off'); 
            else
                set(self.handles.indexEdit, 'enable', 'on'); 
                set(self.handles.slider, 'enable', 'on'); 
                % if (self.value == self.max) | isnan(self.value)
                %     set(self.handles.nextButton, 'enable', 'off');
                % end
                % if (self.value < self.max)
                %     set(self.handles.nextButton, 'enable', 'on');
                % end
                % if (self.value == self.min) | isnan(self.value)
                %     set(self.handles.prevButton, 'enable', 'off');
                % end
                % if (self.value > self.min)
                %     set(self.handles.prevButton, 'enable', 'on');
                % end
            end
        end
    end
end
