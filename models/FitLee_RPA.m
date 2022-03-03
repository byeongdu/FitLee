function out = FitLee_RPA(varargin)
% RPA fit for polydisperse polymer blends.
% Input parameters: 
%   Contrast = (b1/v1 - b2/v2)^2
%               bx : the unit scattering length density of sample x
%               vx : the segment volume of each components.
%               bx/vx : the scattering length density
%   Rgx : The radius of gyration (A) of x
%   chi : Chi parameter
%   vx : the segmental volume (cm^3/mol) of x
%   Nwx : weight averaged degree of polymerization of x
%         if monodisperse, set Nwx = 0;
%   Nnx : number averaged degree of polymerization of x
%   phix : volume fraction of x
%
% Byeongdu Lee
% Ref:
% B. Lee at al. Polymer, 51, 24, 12, 5799–5806.
% B. Lee et al. J. Appl. Cryst.;2009; 42; 161-168.
%
% See also RPAblend.m 
FitLee_helpstr = {'RPA fit for polydisperse polymer blends. ' ,...
'Input parameters:  ',...
'  Contrast = (b1/v1 - b2/v2)^2 ',...
'              bx : the unit scattering length density of sample x ',...
'              vx : the segment volume of each components. ',...
'              bx/vx : the scattering length density ',...
'      Here, it is just a constant.',...
'  RgNx : The radius of gyration (A) of x ',...
'       RgN = N*l/sqrt(6)',...
'       Relation between RgN and Rgz is found in ref 2',...
'  chix1000 : Chi parameter x 1000 ',...
'  vx : the segmental volume (cm^3/mol) of x ',...
'  Nwx : weight averaged degree of polymerization of x ',...
'        if monodisperse, set Nwx = 0 ',...
'  Nnx : number averaged degree of polymerization of x ',...
'  phix : volume fraction of x ',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    1. B. Lee at al. Polymer, 2010, 51, 24, 5799–5806. ',...
'    2. B. Lee et al. J. Appl. Cryst., 2009, 42, 161-168. '};

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
    bestP.contrast = 1;
    bestP.RgnD = 100;
    bestP.RgnH = 45;
    bestP.chix1000 = 3E-1;
    bestP.vD = 32.9;% cm^3/mol
    bestP.vH = 32.92;
    bestP.NwD = 3000;
    bestP.NnD = 1500;
    bestP.NwH = 4000;
    bestP.NnH = 1200;
    bestP.phiD = 0.3;
    bestP.background = 0;
    assignin('base', 'bestP', bestP);
    out = bestP;
    return
end

%% fitting codes .......................
if iscell(q);
    if numel(q) > 1
        error('FitLee_RPA.m is for fitting a set of data, for now')
    end
    q = q{1};
end
q = q(:);


contrast = p.contrast;
RgnD = p.RgnD;
RgnH = p.RgnH;
chi = p.chix1000 / 1000;
vD = p.vD;
vH = p.vH;
NwD = p.NwD;
NnD = p.NnD;
NwH = p.NwH;
NnH = p.NnH;
phiD = p.phiD;
background = p.background;


if NwD ~=0 % polydisperse
    UD = (NwD/NnD-1);
else % monodisperse
    UD = 0;
end
if NwH ~=0
    UH = (NwH/NnH-1); % polydisperse
else
    UH = 0;  % monodisperse
end
xD = q.^2*RgnD^2/(1+2*UD);
xH = q.^2*RgnH^2/(1+2*UH);
if UD ~= 0
    PD = 2*((1+UD*xD).^(1/UD)+xD-1)./(1+UD)./xD.^2;
else
    PD = 2*(exp(-xD)+xD-1)./xD.^2;
end
if UH ~= 0
    PH = 2*((1+UH*xH).^(1/UH)+xH-1)./(1+UH)./xH.^2;
else
    PH = 2*(exp(-xH)+xH-1)./xH.^2;
end

%PD = 2./(xD.^2).*(xD -1+(hD./(hD+xD)).^hD);
%PH = 2./(xH.^2).*(xH -1+(hH./(hH+xH)).^hH);

phiH = 1-phiD;
v =1/(phiD/vD + phiH/vH);
S = 1./(1./(vD*phiD*NnD*PD) + 1./(vH*phiH*NnH*PH) - 2*chi/v);
Intensity = contrast*S;

out = background + Intensity;

if isnan(out)
    out = ones(size(out));
end