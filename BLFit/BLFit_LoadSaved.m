function FitLeehandle = BLFit_LoadSaved(handles, fn)
if nargin < 1
    noFitLee = 1;
    handles = [];
else
    noFitLee = 0;
end

if nargin > 1
    filen = fn;
else
    [filename, pathname] = uigetfile('*.txt; *.dat', 'Pick your file');
%    filen = fullfile(pathname, filename);
    if isequal(filename,0) || isequal(pathname,0)
       disp('User pressed cancel')
       return
    else
       disp(['User selected ', fullfile(pathname, filename)])
       filen = fullfile(pathname, filename);
    end
end

fid = fopen(filen);
rl = fgetl(fid);
if ~ischar(rl); fclose(fid); return; end
[~, rm] = strtok(rl, ':');
FitFunc = rm(3:end);

if isempty(handles)
    noFitLee = 1;
end

if noFitLee
    FitLeehandle = FitLee(FitFunc);
    handles = guihandles(FitLeehandle);
end
while 1
    done = strfind(lower(rl), 'fit parameter');
    if done >0, break, end
    rl = fgetl(fid);
    if ~ischar(rl), break, end
end
[p, qfit, LB, UB, pnames] = BLFit_readparameters(handles);

for i=1:numel(p)
    rl = fgetl(fid);
    if strfind(rl, pnames{i})
        dt = sscanf(rl((numel(pnames{i}))+1:end),'%f');
        p(i) = dt(1);
        qfit(i) = dt(2);
        LB(i) = dt(3);
        UB(i) = dt(4);
    end
end
BLFit_setparameters(handles, p, qfit, LB, UB);

rl = fgetl(fid);
done = strfind(lower(rl), 'fit range');
if done > 0
    rl = fgetl(fid);
    [~, L] = strtok(rl, ' ');
    rl = fgetl(fid);
    [~, R] = strtok(rl, ' ');
    set(handles.edit_Lindex, 'string', L);
    set(handles.edit_Rindex, 'string', R);
end
done = 0;
while ~done
    done = strfind(rl, 'data ======================');
    rl = fgetl(fid);
    if ~ischar(rl), break, end
end

data = [];
rl = fgetl(fid);
while 1
    if ~ischar(rl), break, end
    dt = sscanf(rl, '%f');
    data = [data;dt(:)'];
    rl = fgetl(fid);
end
fclose(fid);
FitLeehandle = get(handles.versionPanel, 'parent');
fit = getappdata(FitLeehandle, 'Fit');
fit.simtag = '__sim';
% try
%     fit = evalin('base', 'fit');
% catch
%     fit = [];
%     fit.simtag = '__sim';
% end
if ~isfield(fit, 'figh')
    fit.figh = figure;
else
    if ~ishandle(fit.figh)
        fit.figh = figure;
    end
end
try
    figure(fit.figh);
    clf(fit.figh);
    loglog(data(:,1), data(:,2), 'bo')
    line('xdata', data(:,1), 'ydata', data(:,3), 'color','r', 'tag', fit.simtag, 'parent', findobj(fit.figh, 'type', 'axes'));
    title(filen)
catch
    figure;
    loglog(data(:,1), data(:,2), 'bo')
    line('xdata', data(:,1), 'ydata', data(:,3), 'color','r', 'tag', fit.simtag, 'parent', findobj(fit.figh, 'type', 'axes'));
    title(filen)
end
fit.fit = data;
fit.filen = filen;
fit.NdataSet = 1;
setappdata(FitLeehandle, 'Fit', fit);
assignin('base', 'fit', fit);