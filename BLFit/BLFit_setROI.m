function fit = BLFit_setROI(fit, xd, yd, handles)

if nargin <4
    handles = fit.handles;
end
noXD = 1;
if nargin > 2
    if ~isempty(xd)
        noXD = 0;
    end
end

if noXD
    try
        xd = get(fit.datahandle(1), 'xdata');
    catch
        error('The plot is changed. Click "Pick Figure" and hit enter')
    end
    yd = get(fit.datahandle(1), 'ydata');
    xd=xd(:);
    yd=yd(:);
end

Numpnt = numel(xd);
neednewdots = 0;
if ~isfield(fit, 'dothandles')
    neednewdots = 1;
else
    if ~ishandle(fit.dothandles(1))
        neednewdots = 1;
    else
        if get(fit.dothandles(1), 'parent') ~= get(fit.figh, 'children')
            neednewdots = 1;
        end
    end
end

if neednewdots
    try
        figure(fit.figh)
    catch
        fit.figh = figure;
    end
    hold on;
    roi = [1, numel(xd)];
    if isfield(fit, 'roi')
        if ~isempty(fit.roi)
            roi = fit.roi;
        end
    end
    roioutofrange = roi > numel(xd);
    roi(roioutofrange) = numel(xd);
    Lh = plot(xd(roi(1)), yd(roi(1)), 'ro', 'markerfacecolor', [.49, 1, .63], 'markersize', 10);
    Rh = plot(xd(roi(2)), yd(roi(2)), 'ro', 'markerfacecolor', [.75, 1, .25], 'markersize', 10);
    set(Lh, 'tag', 'Dot');
    set(Rh, 'tag', 'Dot');
    hold off;
    fit.dothandles = [Lh, Rh];
end


Lind = str2double(get(handles.edit_Lindex, 'string'));
Rind = str2double(get(handles.edit_Rindex, 'string'));
if isnan(Lind)
    Lind = 1;
end
if isnan(Rind)
    Rind = 1;
end

if ~isempty(Lind)
    qL = xd(Lind);
    set(handles.edit_Lq, 'string', num2str(qL));
    set(fit.dothandles(1), 'xdata', xd(Lind)); 
    set(fit.dothandles(1), 'ydata', yd(Lind));
end
if ~isempty(Rind)
    try
        qR = xd(Rind);
    catch
        Rind = Numpnt;
        qR = xd(Rind);
    end
    set(handles.edit_Rq, 'string', num2str(qR));
    set(fit.dothandles(2), 'xdata', xd(Rind)); 
    set(fit.dothandles(2), 'ydata', yd(Rind));
end
fit.roi = [Lind, Rind];