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
dt   = cell(1, numel(filen));
hdrs = cell(1, numel(filen));
for i=1:numel(filen)
    filename = fullfile(pathname, filesep, filen{i});
    [hdrs{i}, data] = hdrload(filename);
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

% If the data header mentions xlabel / ylabel, set the standard SAXS axis labels.
hdrtext = lower(reshape(hdrs{end}.', 1, []));
if contains(hdrtext, 'xlabel') || contains(hdrtext, 'ylabel')
    xlabel('q (\AA^{-1})')
    ylabel('I(q) (cm^{-1})')
end
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

