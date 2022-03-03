function BLFit_save(handles, fn)
isFnNeeded = 0;
if nargin < 2
    isFnNeeded = 1;
    fn = [get(gcbf, 'tag'), '.txt'];
else
    if isempty(fn)
        isFnNeeded = 1;
    end
end
if isFnNeeded
    [filename, pathname] = uiputfile({'*.txt', 'Text Files (*.txt)';...
        '*.dat', 'Ascii Data Files (*.dat)';...
        '*.*', 'All Files (*.*)'}, ...
        'Save as');
    if isequal(filename,0) || isequal(pathname,0)
       disp('User pressed cancel')
       return
    else
       disp(['User selected ', fullfile(pathname, filename)])
       fn = fullfile(pathname, filename);
    end
end

figh = handles.FitLee;
fit = getappdata(figh, 'Fit');

%fit = evalin('base', 'fit');

fid = fopen(fn, 'w');
[p, qfit, LB, UB, pnames] = BLFit_readparameters(handles);
Lindex = get(handles.edit_Lindex, 'string');
Rindex = get(handles.edit_Rindex, 'string');
fitFN = handles.FitLee.Name;
% blk = strfind(fitFN, ' ');
% fitFN = fitFN(blk(end)+1:end);
fprintf(fid, '# Fitting Function: %s\n', fitFN);
fprintf(fid, '# Date of fit: %s\n', datestr(now));

[~,b] = system('echo %username%');
fprintf(fid, '# Who did the fit: %s', b);

[~,b] = system('echo %computername%');
fprintf(fid, '# Computer name: %s', b);

if isfield(fit, 'title')
    fprintf(fid, '# Title of data: %s\n', fit.title);
end
fprintf(fid, '# Fit parameter ============= \n');
for i=1:numel(p)
    fprintf(fid, '%s %0.8e %i %0.8e %0.8e\n', pnames{i}, p(i), qfit(i), LB(i), UB(i));
end
fprintf(fid, '# Fit range ============= \n');
fprintf(fid, '%s %s\n', 'Min', Lindex);
fprintf(fid, '%s %s\n', 'Max', Rindex);

fprintf(fid, '# data ====================== \n');
%for i=1:numel(fit.fit(:,1))
%    fprintf(fid, '%0.8e %0.8e %0.8e\n', fit.fit(i,1), fit.fit(i,end), fit.fit(i,2));
%end
dt = [fit.fit(:,1), fit.fit(:,end), fit.fit(:, 2:end-1)];
fclose(fid);
dlmwrite(fn, dt,'delimiter','\t','precision','%.8f', '-append');

