function fit = BLFit_qselect(varargin)
hdl = [];indx = [];xd=[];yd=[];zd=[];tg=[];scanIndx=[];
FitLeehandle = [];
if numel(varargin) > 0
    if ishandle(varargin{1})
        FitLeehandle = varargin{1};
    end
end
if ~isempty(FitLeehandle)
    fit = getappdata(varargin{1}, 'Fit');
else
    fit = evalin('base', 'fit');
end

Figh = fit.figh;
ax = findobj(Figh, 'type', 'axes');
t = findobj(Figh, 'Tag', 'Dot');
if (numel(t) > 2)
    if strcmp(get(ax, 'yscale'), 'log');
        k = cell2mat((get(t, 'ydata')))==0;
        t(k) = [];
    else
        error('There are more than two points are selected')
    end
end

if numel(t) ~=2
    error('Should select only two points ')
end
    
        
if isempty(t)
    disp('error!! you did not select any curve');
    return
end    

dothdl = t;

%if ~ishandle(fit.datahandle)
%    fit = [];
%    fit.figh = Figh;
%    fit.datahandle = findobj(Figh, 'Tag', '', 'type', 'line');
%    fit.xd = get(fit.datahandle, 'xdata');
%    fit.yd = get(fit.datahandle, 'ydata');
%    fit.simlineh = findobj(Figh, 'Tag', 'simulation', 'type', 'line');
%    fit.fitlineh = [];
%end

for i=1:numel(dothdl)
    tmpd = findobj(get(dothdl(i), 'userdata'), 'type', 'line');
    Ndata = get(tmpd, 'tag');
    Ndata = str2num(Ndata(isnumber(Ndata)));
    %dthandle(i) = findobj(Figh, 'tag', fit.dattag, 'userdata', get(tmpd, 'userdata'));
    dthandle(i) = fit.datahandle(Ndata);
end
%dthandle = fit.datahandle;
[dthdle, ~, dthidx] = unique(dthandle);
if numel(dothdl) ~= 2*numel(dthdle)
    error('Selection wrong, see BLFit_qselect.m')
end

roi = [];
data2fit = [];
numpnts = [];
for i=1:numel(dthdle)
    dot = dothdl(dthidx == i);
    if numel(dot) ~= 2
        error('there should be two dots per data')
    end
    for m=1:numel(dot)
        indx(m) = getindx(dot(m), dthdle(i));
    end
    indx = sort(indx);
    roi = [roi; indx];
    xd = get(dthdle(i), 'xdata');xd=xd(:);
    yd = get(dthdle(i), 'ydata');yd=yd(:);
    data = [xd, yd];
    fit.data{i} = data;
    xdt = xd(indx(1):indx(2));%xdt = xdt(:);
    ydt = yd(indx(1):indx(2));%ydt = ydt(:);
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
    lh = line('xdata', data2fit(sidx:eidx,1), 'ydata', data2fit(sidx:eidx,2), 'color','r', 'parent', findobj(Figh, 'type', 'axes'));
    fitlh = [fitlh, lh];
end
set(fitlh, 'tag', 'fit');

fit.data2fit = data2fit;
fit.numpnts = numpnts;
fit.roi = roi;
fit.fitlineh = fitlh;
if ~isempty(FitLeehandle)
    setappdata(FitLeehandle, 'Fit', fit);
end

assignin('base', 'fit', fit);

function [indx, xr, yr, xv, yv] = getindx(t, handle)

    xd = get(t, 'Xdata');
    yd = get(t, 'Ydata');
    xv = get(handle(1), 'xdata');
    yv = get(handle(1), 'ydata');
    indx = zeros(2,1);
    xr=zeros(2,1);yr=zeros(2,1);
    if numel(xd) == 1
        indx = find(xv == xd);
        xr = xv(indx);
        yr = yv(indx);
        return
    end

    indx = find(xv == xd);
    xr = xv(indx);
    yr = yv(indx);
