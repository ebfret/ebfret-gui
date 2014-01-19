function [series, analysis] = select_analysis(parent)
    function resume(int)
        status = int;
        uiresume();
    end

    % function update(par, field, value)
    %     if nargin < 1
    %         pars = fieldnames(values);
    %         for p = 1:length(pars)
    %             fields = fieldnames(self.edit.(pars{p}));
    %             for f = 1:length(fields)
    %                 update(pars{p}, fields{f}, ...
    %                     num(pars{p}, fields{f})); 
    %             end
    %         end
    %     else
    %         values.(par).(field) = value;
    %     end
    % end

    % function refresh()
    %     pars = fieldnames(self.edit);
    %     for p = 1:length(pars)
    %         fields = fieldnames(self.edit.(pars{p}));
    %         for f = 1:length(fields)
    %             try
    %                 val = ebfret.num_to_str(values.(pars{p}).(fields{f}));
    %                 set(self.edit.(pars{p}).(fields{f}), 'string', val{1});
    %             catch
    %             end
    %         end
    %     end
    % end 

    function value = num(par, field)
        value = str2num(get(self.edit.(par).(field), 'string'));
    end

    % global vars
    status = 0;
    values = struct();
    self = struct();

    self.dialog = dialog('name', 'Select', ...
                         'color', [0.95 0.95 0.95], ...
                         'units', 'pixels', ...
                         'CloseRequestFcn', @(varargin) resume(0));
    set(self.dialog, ...
        'DefaultUIPanelBackGroundColor', [0.95 0.95 0.95], ...
        'DefaultUIControlUnits', 'normalized');

    try
        % element sizes (pixels)
        row_height = 18;
        edit_width = 72;
        label_width = 72;

        % vertical and horizontal padding (pixels)
        pad_width = 6;
        pad_height = 6;

        % dialog height and width
        num_rows = 5;
        dialog_height = (num_rows+1) * row_height + (num_rows+2) * pad_height;
        dialog_width = 3 * pad_width + label_width + edit_width;

        % get normalized units
        rh = row_height / dialog_height;
        ew = edit_width / dialog_width;
        lw = label_width / dialog_width;

        ph = pad_height / dialog_height;
        pw = pad_width  /dialog_width;

        % adjust dialog size
        rect = get(self.dialog, 'position');
        set(self.dialog, 'position', [rect(1) rect(2) dialog_width dialog_height]);

        % get popup menu options
        states = num2cell(parent.controls.min_states:parent.controls.max_states);
        groups = cat(2, {'all'}, sort(unique({parent.series.group})));

        self.label.analysis ...
            = uicontrol(self.dialog, ...
                'style', 'text', ...
                'backgroundcolor', [0.95 0.95 0.95], ...
                'string', 'Number of States', ...
                'horizontalalignment', 'center', ...
                'position', [pw 1-1*(rh+ph) 1-2*pw rh]);

        self.popup.analysis ...
            = uicontrol(self.dialog, ...
                'style', 'popup', ...
                'string', ebfret.join('|', cellfun(@num2str, states, 'UniformOutput', false)), ...
                'position', [pw 1-2*(rh+ph) 1-2*ph rh]);

        self.label.group ...
            = uicontrol(self.dialog, ...
                'style', 'text', ...
                'backgroundcolor', [0.95 0.95 0.95], ...
                'string', 'Time Series Group', ...
                'horizontalalignment', 'center', ...
                'position', [pw 1-3*(rh+ph) 1-2*pw rh]);

        self.popup.group ...
            = uicontrol(self.dialog, ...
                'style', 'popup', ...
                'string', groups, ...
                'position', [pw 1-4*(rh+ph) 1-2*ph rh]);

        self.okButton ...
            = uicontrol(self.dialog, ...
                'style', 'pushbutton', ...
                'string', 'Ok', ...
                'position', [pw ph 0.5-1.5*pw rh], ...
                'callback', @(varargin) resume(1));
                % 'callback', ...
                %     @(source, event) return_values(1), ...
        self.cancelButton ...
            = uicontrol(self.dialog, ...
                'style', 'pushbutton', ...
                'string', 'Cancel', ...
                'position', [0.5+0.5*pw ph 0.5-1.5*pw rh], ...
                'callback', @(varargin) resume(0));

        %methodPopupCallback(0);
        % refresh();
        uiwait(self.dialog);
        state = states{get(self.popup.analysis, 'value')};
        analysis = parent.analysis(state);
        group = groups{get(self.popup.group, 'value')};
        if strcmpi(group, 'all')
            series = parent.series;
        else
            ns = find(arrayfun(@(s) strcmpi(s.group, group), parent.series));
            series = parent.series(ns);
            analysis.viterbi = analysis.viterbi(ns);
            analysis.expect = analysis.expect(ns);
            analysis.posterior = analysis.posterior(ns);
            analysis.lowerbound = analysis.lowerbound(ns);
            analysis.restart = analysis.restart(ns);
            analysis.prior = ebfret.analysis.hmm.h_step(analysis.posterior, analysis.prior);
        end
        delete(self.dialog);
    catch err
        delete(self.dialog);
        rethrow(err)
    end
end