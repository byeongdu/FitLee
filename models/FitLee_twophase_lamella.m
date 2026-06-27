function [out, report] = FitLee_twophase_lamella(varargin)
    % Paracrystal theory equation for triblock copolymer.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Two-phase 1D lamella stack diffraction.' ,...
'',...
'$I(q) = I_0 \cdot $ lam_twophase$(q, [d, v_f], \Delta\rho, domainsize, microstrain, DW) + I_b$',...
'$\qquad    I_b = poly1 \cdot q + poly2$',...
'',...
'Periodic stack with repeat distance $d$ and phase-A thickness $v_f \cdot d$.',...
'Peak positions at $q_n = 2\pi n / d$.',...
'',...
'$Parameters$',...
'$\quad  I0$ : overall scale factor',...
'$\quad  d$ : lamella repeat distance (A)',...
'$\quad  vf$ : volume fraction (== thickness fraction) of phase A',...
'$\quad  drho$ : electron-density contrast between A and the other phase',...
'$\quad  domainsize$ : stack coherence length (A)',...
'$\quad  microstrain$ : RMS strain broadening of the layer spacing',...
'$\quad  DW$ : Debye-Waller factor',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    Senesi, A. J. and B. Lee. Small-Angle Scattering of Particle Assemblies.',...
'    J. Appl. Cryst. (2015) 48(4): 1172-1182.'};

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
    bestP.I0 = 1;
    bestP.d = 200;
    bestP.vf = 0.3;
    bestP.drho = 0.1;
    bestP.domainsize = 5000;
    bestP.microstrain = 0.02;
    bestP.DW = 0.1;
    bestP.poly1 = 0;
    bestP.poly2 = 0;
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
q = q(:);

Iq = lam_twophase(q, [p.d, p.vf], p.drho,p.domainsize,p.microstrain, p.DW);

back = p.poly1*q + p.poly2;
out = p.I0*Iq + back;

if isnan(out)
    out = ones(size(out))*1E20;
end

end
