function [ret, rtn] = findcellstr(cstr, fstring)
% [ret, rtn] = findcellstr(cstr, fstring)
% findcellstr is to find the string, fstring, in the cell, cstr.
% and return index or indice..
% rtn is 1 or 0
%
% cstr : cellstring
% fstring : string to find
% Byeongdu Lee, Nov 28, 2007

if ischar(fstring)
    [ret, rtn] = getindex(cstr, fstring);
elseif iscell(fstring)
    ret = zeros(size(fstring));
    for i=1:numel(fstring)
        ret(i) = getindex(cstr, fstring(i));
    end
end

function [ret, rtn] = getindex(cstr, fstr)
rtn = strcmp(cstr, fstr);
ret = find(rtn == 1);
return

% old code ...
ret = [];
for i=1:numel(cstr)
    if strcmp(cstr{i}, fstring)
        ret = [ret, i];
    end
end
if isempty(ret)
    ret = 0;
end