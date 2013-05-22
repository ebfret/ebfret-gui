function str = num_to_str(num, max_length)
    if nargin < 2
        max_length = 4;
    end
    function str = strip_zero(str)
        ind = find(str ~= '0', 1, 'last');
        if str(ind) == '.'
            ind = ind - 1;
        end
        str = str(1:ind);
    end
    function str = format_sci(num)
        if num == 0
            str = '0';
            return
        end
        dec = floor(log(num) ./ log(10));
        fl = (num ./ 10^dec);
        if dec ~= 0
            if round(fl) == fl
                str = sprintf('%de%d', fl, dec);
            else
                str = sprintf('%.1fe%d', fl, dec);
            end
        else
            if round(fl) == fl
                str = sprintf('%d', fl);
            else
                str = sprintf('%.1f', fl);
            end
        end
    end
    % try floating point notation first
    str = arrayfun(@(num) strip_zero(sprintf('%f', num)), num, 'UniformOutput', false);
    if any(cellfun(@length, {str{:}}) > max_length)
        % return scientific notation if float notation is too long
        str = arrayfun(@(num) format_sci(num), num, 'UniformOutput', false);
    end
end
