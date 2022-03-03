function fit = BLFit_index_callback(varargin)
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

%handles = fit.handles;

fit = BLFit_setROI(fit, [], [], handles);
if ~isempty(FitLeehandle)
    setappdata(FitLeehandle, 'Fit', fit);
end
assignin('base', 'fit', fit);