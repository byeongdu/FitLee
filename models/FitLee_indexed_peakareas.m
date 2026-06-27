function [out, report] = FitLee_indexed_peakareas(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
global FitLee_helpstr
FitLee_helpstr = {'Fit peak areas for a curve indexed by indexing ' ,...
'',...
'Note that',...
'I(q)_diff_cross_section = (delta\_rho*r_e)^2*f_n*P(q)',...
'         , where P(0) = V_p^2'...
' ',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    1. Kwon et al. Nature Materials, 2015, 14, 215–223. ',...
'    2. Wang et al. J. Phys. Chem. C. 2013. 117(44), 22627. '};

if numel(varargin) > 1
    p = varargin{1};
    q = varargin{2};
    isini = 0;
elseif isscalar(varargin)
    p = varargin{1};
    isini = 1;
    if ishandle(p)
        guiadd(p);
        return
    end
elseif numel(varargin) == 0
    peaks = evalin('base', 'indexing_hkls');
    p = size(peaks,1);
    isini = 1;
end

%% initialize fit parameter bestP
if isini
    if exist('pseudovoigt', 'file') ~= 2
        error('Need "pseudovoigt.m". Download and keep in with FitLee_multipseudovoigt.m.')
    end
    Nf = p;
    bestP = [];

    for i=1:Nf
        bestP.(['Area', num2str(i)]) = peaks(i, 2);
    end
    % Need 4 parameters for background.
    bestP.Gaussian_Width = 0.005;
    bestP.Domain_Size = 1000;
    bestP.Microstrain = 0.001;
    bestP.SF_userBG = 1;
    bestP.string = 'FitLee_indexed_peakareas(figH)';
    assignin('base', 'bestP', bestP);
    out = bestP;
    return
end

%% fitting codes .......................
if iscell(q)
    if numel(q) > 1
        error('FitLee_schultzsphere.m is for fitting a set of data, for now')
    end
    q = q{1};
end

back_goodpix = [];

try
    UBG = evalin('base', 'userbackground');
    q = round(q, 6);
    UBG(:,1) = round(UBG(:,1), 6);
    k = UBG(:,1)<q(1) | UBG(:,1)>q(end);
    UBG(k,:) = [];
    if numel(q) == numel(UBG(:,1))
        if size(UBG, 2)==3
            back_goodpix = logical(UBG(:,3));
        else
            back_goodpix = logical(size(UBG(:,2)));
        end
        bg = UBG(:,2);
    else
        bg = interp1(UBG(:,1), UBG(:,2), q);
    end
catch
    bg = zeros(size(q));
end

if ~isempty(back_goodpix)
    dt_goodpix = ~back_goodpix;
    p_start = find(diff(dt_goodpix)==1)+1;
    p_end = find(diff(dt_goodpix)==-1);
else
    p_start = [];
    p_end = [];
end
peaks = evalin('base', 'indexing_hkls');
back = p.SF_userBG*bg;

predYs = {};
predY = zeros(size(q));
pnames = fieldnames(p);
np = 0;
for i=1:numel(pnames)
    if strfind(pnames{i}, 'Area')
        np = np + 1;
    end
end

for k=1:np
    w =  peakwidth(peaks(k,1), 0.5, p.Domain_Size, p.Microstrain);
    qv = peaks(k,1);
    pa = [p.(['Area', num2str(k)]), ...
        qv, p.Gaussian_Width, w, 0];
    pa = abs(pa);
    y = abs(pseudovoigt(q, pa));
    if ~isempty(p_start)
        [~, ind] = min(abs(q-qv));
        t = find(p_start<ind);pind = t(end);
        p2 = p_end(pind);
        p1 = p_start(pind);
        pn = polyfit([q(p1), q(p2)], [y(p1), y(p2)], 1);
        y = y-(pn(1)*q+pn(2));
        y(1:(p1-1)) = 0;
        y((p2+1):end) = 0;
    end
    predY = predY + y;
end
predY = predY + back;
out = predY;
if isnan(out)
    out = ones(size(out));
end

function guiadd(figH)
    uph = findobj(figH, 'type', 'uipanel');
    pos = get(figH, 'position');
    hFigHeight = pos(end);
        uicontrol(...
          'style'                     , 'text', ...
          'backgroundcolor'           , [0,0.5,0.5], ...
          'foregroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'Left', ...
          'parent'                    , uph, ...
          'string'                  ,   'BackG',...
          'position'                  ,[1, hFigHeight-210,100,20], ...
          'tag'                       , 'text_Pshape');
    popup = uicontrol(...
          'parent', uph,...
          'horizontalalignment'       , 'Left', ...
          'Style', 'pushbutton',...
          'String', 'Load Background',...
          'Position', [50, hFigHeight-235, 150, 50],...
          'Callback', @setshape);
function setshape(varargin)
    obj = varargin{1};
    [filename, pathname] = uigetfile('*.txt; *.dat; *.*', 'Pick your file');
    if isequal(filename,0) || isequal(pathname,0)
        disp('User pressed cancel')
        return
    else
        disp(['User selected ', fullfile(pathname, filename)])
    end
    filename = fullfile(pathname, filesep, filename);
    [~, data] = hdrload(filename);
    assignin('base', 'userbackground', data)