function [out, report] = FitLee_powerlaw(varargin)
FitLee_helpstr = {'Powerlaw fit. ' ,...
' y = a*(x+xc)^alpha',...
'   a : amplitude',...
'   xc : position',...
'   alpha : exponent', ...
'',...
'Byeongdu Lee (blee@anl.gov)'};

if numel(varargin) > 1
    p = varargin{1};
    q = varargin{2};
    isini = 0;
elseif numel(varargin) == 1
    p = varargin{1};
    isini = 1;
    if ischar(p)
        out = FitLee_helpstr;
        return
    end
elseif numel(varargin) == 0
    p = 1;
    isini = 1;
end

%% initialize fit parameter bestP
if isini
    bestP = [];
    bestP.a = 5;
    bestP.x0 = 273.15;
    bestP.alpha = -1/2;
    assignin('base', 'bestP', bestP);
    out = bestP;
    return
end

%% fitting codes .......................
if iscell(q)
    if numel(q) > 1
        error('This function is for fitting a set of data, for now')
    end
    q = q{1};
end

q = q(:);
y = p.a*(q+p.x0).^(p.alpha);
out = y;
if isnan(out)
    out = ones(size(out));
end

if nargout == 2
    %q = linspace(1E-5, 1, 2);
    report = '';
end

