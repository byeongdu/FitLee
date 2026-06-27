function out = FitLee_rocking(varargin)
FitLee_helpstr = {'Rocking-scan V-shape fit (piecewise sine of angle).',...
'',...
'$y(\theta) = \left\{ \begin{array}{ll} +L_u \cdot \sin((\theta - a_0)\pi/180) + z_p & \theta > a_0 \\ -L_d \cdot \sin((\theta - a_0)\pi/180) + z_p & \theta \leq a_0 \end{array} \right.$',...
'',...
'$\theta$ is in degrees. Curve forms a V (or inverted V) with vertex at $(a_0, z_p)$',...
'and slopes $L_u$ on the right and $L_d$ on the left.',...
'',...
'$Parameters$',...
'$\quad  a0$ : vertex angle (degree)',...
'$\quad  Lu$ : upper-branch amplitude',...
'$\quad  Ld$ : lower-branch amplitude',...
'$\quad  zp$ : vertical offset of the vertex',...
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
