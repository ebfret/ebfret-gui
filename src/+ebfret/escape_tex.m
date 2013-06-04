function str = escape_tex(str)
    % escapes TeX characters in a string to prevent matlab rendering.
    str = regexprep(str, '([\^_{}]{1})', '\\$1');