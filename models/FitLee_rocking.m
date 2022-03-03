function out = FitLee_rocking(varargin)
FitLee_helpstr = {'\bf{Rocking scan fit}',...
'\rm This function is only for fitting a set of data at a time.',...
'if numel(varargin) == 1 and p is not a struct but a number, then',...
'it generate set of default parameters for using this function.',...
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
    %Nf = p;
    bestP = [];
    bestP.a0 = 0.001;
    bestP.Lu = 2.5;
    bestP.Ld = 2.5;
    bestP.zp = -0.001;
%    bestP.zp = 0;
    out = bestP;
    return
end

if iscell(q);
    if numel(q) > 1
        error('FitLee_rocking.m is for fitting a set of data, for now')
    end
    q = q{1};
end
q = q(:);
t = q > p.a0;
out = zeros(size(q));
if isfield(p, 'zp');
    zp = p.zp;
else
    zp = 0;
end
out(t) = p.Lu*sin((q(t)-p.a0)*pi/180)+zp;
t = q <= p.a0;
out(t) = -p.Ld*sin((q(t)-p.a0)*pi/180)+zp;
