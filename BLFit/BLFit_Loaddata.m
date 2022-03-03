function BLFit_Loaddata(handles)
%try
%    fit = evalin('base', 'fit');
%catch
%    fit = [];
%end
[filename, pathname] = uigetfile('*.txt; *.dat', 'Pick your file', 'MultiSelect', 'on');
n = figure;
if ~iscell(filename)
    filen{1} = filename;
else
    filen = filename;
end
cl = ['o','^','v','s','+'];
dt = {};
for i=1:numel(filen)
    filename = fullfile(pathname, filesep, filen{i});
    [~, data] = hdrload(filename);
    k = isnan(data(:,1));data(k,:) = [];
    k = isnan(data(:,2));data(k,:) = [];
    dt{i} = data;
end
%if size(dt{1}, 2) > 2
if ishandle(handles.FitLee)
    pstyle = getappdata(handles.FitLee, 'plotstyle');
end
if ~isempty(pstyle)
    Xcol = pstyle.Xcol;
    Ycol = pstyle.Ycol;
    pstyle = pstyle.plotstyle;
else
    Xcol = 1;
    Ycol = 2;
    pstyle = 'plot';
end
for i=1:numel(filen)
    t = feval(pstyle, dt{i}(:,Xcol), dt{i}(:,Ycol), ['b',cl(i)]);hold on;
%    t = loglog(dt{i}(:,1), dt{i}(:,2), ['b',cl(i)]);hold on;
%    fit.datahandle(i) = t;
    filename = fullfile(pathname, filesep, filen{i});
    set(t, 'tag', filename);
end
titlestr = filen{end};
t = strfind(titlestr, '_');titlestr(t) = ' ';
title(titlestr);
warning off
matlabver = version;
t = strfind(matlabver, '.');
mver = str2double(matlabver(1:t(2)-1));
if mver <=8.3
    set(handles.edit_Figure, 'string', num2str(n));
else
    set(handles.edit_Figure, 'string', num2str(n.Number));
end
%assignin('base', 'fit', fit)
end

