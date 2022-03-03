function BLFit_pickFigure(Figh, handles)
% Function to determine the fit data.
% Both single data set and multiple data sets are acceptable.
% fit, a struct variable, will be saved to the base workspace.
% If multiple sets are selected, xd, yd and data fields in fit will be a
% cell array
% field of fit.
%   handles
%   figh
%   NdataSet : the number of fit data
%   xd : it will be a cell array when NdataSet > 1
%   yd : it will be a cell array when NdataSet > 1
%   data : A cell array no matter when NdataSet > 1 or not.
   
% B. Lee 2/8/2014
try
    fit = evalin('base', 'fit');
catch
    fit = [];
    fit.dattag = '__data';
    fit.fittag = '__fit';
    fit.simtag = '__sim';
end
%fit.handles = handles;
fit.figh = Figh;
fit.data = [];
fit.data2fit = [];
fit.tag = {};
fit.datahandle = [];

axish = findobj(Figh, 'type', 'axes');
delete(findobj(axish, 'tag', 'temporary'));
if isfield(fit, 'fittag')
    delete(findobj(axish, '-regexp', 'tag', fit.fittag));
end
if isfield(fit, 'simtag')
    delete(findobj(axish, '-regexp', 'tag', fit.simtag));
end
delete(findobj(axish, 'tag', 'Dot'));

axistitle = get(axish, 'title');
axistitlestr = get(axistitle, 'string');
xd = [];
yd = [];m=1;
t = findobj(get(axish, 'children'), 'type', 'line');
if numel(t) > 1
    rmdata = [];
    for i=1:numel(t)
        if isequal(get(t, 'color'), [0, 1, 0]) | isequal(get(t, 'color'), [1, 0, 0])
            rmdata = [rmdata, i];
        end
    end
    t(rmdata) = [];
end

if isempty(t)
    hFigSAXSLee = findall(0,'Tag','SAXSLee_Fig');
    hAxes = findall(hFigSAXSLee,'Tag','SAXSLee_Axes');
    t = findobj(hAxes, 'type', 'line', '-not', {'-regexp', 'tag', 'BACK'});
end

%added 3/21/2014
if ~isfield(fit, 'fittag')
    fit.fittag = '';
end
if ~isfield(fit, 'simtag')
    fit.simtag = '';
end
if ~isfield(fit, 'dattag')
    fit.dattag = '';
end
% 

for i=1:numel(t);
    tagn = get(t(i), 'tag');
    isfit = ~isempty(strfind(tagn, fit.fittag));
    issim = ~isempty(strfind(tagn, fit.simtag));
    %issim = isempty(strfind(tagn, fit.simtag));
    if ~isfit && ~issim;
        fit.datahandle(m) = t(i);
        if numel(get(t(i), 'xdata')) < 2
            continue
        end
        xd{m} = get(t(i), 'xdata');
        yd{m} = get(t(i), 'ydata');

        data = [xd{m}(:), yd{m}(:)];
        fit.data{m} = data;
        fit.tag{m} = get(t(i), 'tag');
        set(t(i), 'color', 'b');
        set(t(i), 'tag', [fit.dattag, num2str(i)]);
        m = m+1;
    else
        delete(t(i));
    end
end
%cla(axish);
%figure(Figh);
if numel(xd)>1
    if iscell(xd)
        [indx,tf] = listdlg('PromptString',{'Select a data.',''},...
    'SelectionMode','Multiple','ListString',fit.tag);
        if isempty(indx)
            return
        end
        xd = xd(indx);
        yd = yd(indx);
        fit.data = fit.data(indx);
        fit.tag = fit.tag(indx);
        fit.datahandle = fit.datahandle(indx);
    end
%     rep = ones(size(xd));
%     rep(indx) = 0;

%     rep = zeros(size(xd));
%     for i=1:numel(xd)
%         if ~strcmp(fit.tag{i}, 'Dot')
%             reply = input(sprintf('The size of the data %s is %i, Do you like to fit N[Y/N]:', fit.tag{i}, numel(xd{i})),'s');
%         else
%             reply = 'N';
%         end
%         
%         if ~isempty(reply)
%             if strcmp(reply, 'Y')
%                 reply = 0;
%             else
%                 reply = 1;
%             end
%         else
%             reply = 1;
%         end
%         rep(i) = reply;
%     end
%     rep = find(rep == 1);
%     xd(rep) = [];
%     yd(rep) = [];
%     fit.data(rep) = [];
%     fit.tag(rep) = [];
%     fit.datahandle(rep) = [];
end

fit.NdataSet = numel(xd);

if numel(xd) == 1
    xd = xd{1};
    yd = yd{1};
    fit.tag = fit.tag{1};
    %fit.datahandle = loglog(xd, yd);
    %set(fit.datahandle, 'tag', fit.dattag);
    %set(fit.datahandle, 'userdata', 1);
    %simlh = line('xdata', xd, 'ydata', yd, 'color','g', 'parent', axish);
    %set(simlh, 'tag', fit.simtag);
    %set(simlh, 'userdata', 1);
    %fit.simlineh = simlh;
else
    fit.datahandle = datlh;
    %fit.simlineh = simlh;
end
fit.figh = Figh;
% if strcmpi(get(Figh, 'tag'), 'Load_FitLee')
%     fit.figh = figure;
%     fit.datahandle = plot(xd, yd, 'bo-');
%     set(fit.datahandle, 'tag', fit.dattag);
% else
%     fit.figh = Figh;
% end
%title(axistitlestr);

%datalabel

fit.xd = xd;
fit.yd = yd;
%fit.datahandle = lh;
fit.name = get(gcbf, 'tag');
fit = BLFit_setROI(fit, xd, yd, handles);
if nargin > 1
    fit.handles = handles;
end
assignin('base', 'fit', fit)
try
    setappdata(fit.handles.FitLee, 'Fit', fit)
catch
end