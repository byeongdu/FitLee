function out = FitLee_unifiedSAXS(varargin)

FitLee_helpstr = {'Unified fit. ' ,...
'Input parameters:  ',...
'  2 Level of unified eqns',...
'  G, Rg, B, P(>0) ',...
'  y = G*exp(-x.^2*Rg^2/3) + B*(erf(x*Rg/sqrt(6)).^3./x).^P',...
'  G and Rg are I(0) and Rg for Guinier function',...
'  B and P are I(0) and Porod slope for Porod function',... 
'  I(q) = y1 + y2',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    Unified model : G. Beaucage, J. appl. Cryst(1995), 28, 717-728. '};

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
    Nf = p;
    bestP = [];
    bestP.G1 = 100;
    bestP.Rg1 = 500;
    bestP.B1 = 0.0005;
    bestP.P1 = 2;
    bestP.G2 = 100;% cm^3/mol
    bestP.Rg2 = 50;
    bestP.B2 = 0.0001;
    bestP.P2 = 2;
    bestP.background = 0;
    bestP.SF_userBG = 1;
    assignin('base', 'bestP', bestP);
    out = bestP;
    return
end

%% fitting codes .......................
if iscell(q);
    if numel(q) > 1
        error('FitLee_unifiedSAXS.m is for fitting a set of data, for now')
    end
    q = q{1};
end
q = q(:);


try
    UBG = evalin('base', 'userbackground');
    UBG = interp1(UBG(:,1), UBG(:,2), q);
catch
    UBG = zeros(size(q));
end

G1 = p.G1;
Rg1 = p.Rg1;
B1 = p.B1;
P1 = p.P1;
G2 = p.G2;
Rg2 = p.Rg2;
B2 = p.B2;
P2 = p.P2;
    p.I0 = G1;
    p.B = B1;
    p.Rg = Rg1;
    p.P = P1;
    p.D = 0;
    p.vf = 0.0;
    p.powI = 0;
    p.PorodExp = -4;
    % Need 4 parameters for background.
    p.poly1 = 0;
    p.poly2 = -2;
    p.poly3 = 0;
    p.poly4 = 0;

Iq1 = unifiedSAXS(p, q);

    p.I0 = G2;
    p.B = B2;
    p.Rg = Rg2;
    p.P = P2;
    p.D = 0;
    p.vf = 0.0;
    p.powI = 0;
    p.PorodExp = -4;
    % Need 4 parameters for background.
    p.poly1 = 0;
    p.poly2 = -2;
    p.poly3 = 0;
    p.poly4 = 0;

Iq2 = unifiedSAXS(p, q);

background = p.background;
out = Iq1 + Iq2 + background + p.SF_userBG*UBG;

if isnan(out)
    out = ones(size(out));
end