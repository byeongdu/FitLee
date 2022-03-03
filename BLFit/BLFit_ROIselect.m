function BLFit_ROIselect
hdl = [];indx = [];xd=[];yd=[];zd=[];tg=[];scanIndx=[];
fit = evalin('base', 'fit');
Figh = fit.figh;
handles = fit.handles;
roi = [];
data2fit = [];
numpnts = [];
Lind = str2double(get(handles.edit_Lindex, 'string'));
Rind = str2double(get(handles.edit_Rindex, 'string'));
if Lind >= Rind
    error('Left index should be smaller than Right index')
end
roi = [Lind, Rind];

for i=1:numel(fit.datahandle)
    xd = get(fit.datahandle(i), 'xdata');xd=xd(:);
    yd = get(fit.datahandle(i), 'ydata');yd=yd(:);
    data = [xd, yd];
    fit.data{i} = data;
    xdt = xd(Lind:Rind);%xdt = xdt(:);
    ydt = yd(Lind:Rind);%ydt = ydt(:);
    data2fit = [data2fit; xdt, ydt];
    numpnts = [numpnts, numel(xdt)];
end

%for i=1:numel(t);
%    [indx, xd, yd, xv, yv] = getindx(dothdl(i), fit.datahandle);
%    data(:,1) = xv(:);
%    data(:,2) = yv(:);
%    if ~isempty(zd);
%        data(:,3) = zd(:);
%    end
%end
%indx = sort(indx);
%fit.roi = indx;
if isfield(fit, 'fitlineh')
    try
        delete(fit.fitlineh);
    catch
        fit.fitlineh = [];
    end
end

fitlh = [];
eidx = 0;
for m=1:numel(numpnts)
    if m==1
        sidx = 1;
    else
        sidx = eidx + 1;
    end
    eidx = eidx+numpnts(m);
    lh = line('xdata', data2fit(sidx:eidx,1), 'ydata', data2fit(sidx:eidx,2), ...
        'tag', fit.fittag, 'color','r', 'parent', findobj(Figh, 'type', 'axes'));
    fitlh = [fitlh, lh];
end
set(fitlh, 'tag', 'fit');

fit.data2fit = data2fit;
fit.numpnts = numpnts;
fit.roi = roi;
fit.fitlineh = fitlh;

assignin('base', 'fit', fit);