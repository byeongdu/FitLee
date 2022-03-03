function fit = BLFit_qindex_callback(varargin)

FitLeehandle = [];
if numel(varargin) > 0
    if ishandle(varargin{1})
        FitLeehandle = varargin{1};
    end
end
if ~isempty(FitLeehandle)
    fit = getappdata(FitLeehandle, 'Fit');
    handles = guihandles(FitLeehandle);
else
    fit = evalin('base', 'fit');
    handles = fit.handles;
end
if isempty(fit)
    return
end


Lq = str2double(get(handles.edit_Lq, 'string'));
Rq = str2double(get(handles.edit_Rq, 'string'));
if isnan(Lq)
    Lq = 1;
end
if isnan(Rq)
    Rq = 1;
end

try
    xd = get(fit.datahandle(1), 'xdata');
catch
    error('The plot is changed. Click "Pick Figure" and hit enter')
end
yd = get(fit.datahandle(1), 'ydata');
xd=xd(:);
yd=yd(:);

Lind = find(abs(xd-Lq)==min(abs(xd-Lq)));
Rind = find(abs(xd-Rq)==min(abs(xd-Rq)));

if numel(Lind)>1
    Lind = Lind(1);
end
if numel(Rind)>1
    Rind = Rind(1);
end
set(handles.edit_Lindex, 'string', Lind)
set(handles.edit_Rindex, 'string', Rind)

fit = BLFit_setROI(fit, xd, yd);
if ~isempty(FitLeehandle)
    setappdata(FitLeehandle, 'Fit', fit)
end
% if ishandle(varargin{1})
%     setappdata(varargin{1}, 'Fit', fit);
% end

assignin('base', 'fit', fit);