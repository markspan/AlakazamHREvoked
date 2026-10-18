function code = normaliseCode(type)
%NORMALISECODE  An event code as two recordings can share it: spaces gone,
%   a leading S before digits dropped (BrainVision's "S  1" is 1, as
%   DefineBins reads it; a response's R is kept), and a number written
%   without leading zeros ("007" is 7).
%
%   See also SYNC.MATCHTRIGGERS.
    code = regexprep(char(string(type)), '\s', '');
    if ~isempty(regexp(code, '^[Ss]\d+$', 'once'))
        code = code(2:end);
    end
    if ~isempty(regexp(code, '^\d+$', 'once'))
        code = sprintf('%d', str2double(code));
    end
end
