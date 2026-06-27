function [out, report] = FitLee_twophase_HEX(varargin)
    % Paracrystal theory equation for triblock copolymer.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Two-phase 2D hexagonal lattice diffraction (cylinders in matrix).' ,...
'',...
'$I(q) = I_0 \cdot $ hex_twophase$(q, [d, v_f], \Delta\rho, domainsize, microstrain, DW) + I_b$',...
'$\qquad    I_b = poly1 \cdot q + poly2$',...
'',...
'2D hexagonal lattice with edge $d$; unit-cell area $A = d^2 \sqrt{3}/2$.',...
'Implied cylinder radius $R = \sqrt{v_f \cdot A / \pi}$.',...
'',...
'$Parameters$',...
'$\quad  I0$ : overall scale factor',...
'$\quad  d$ : lattice constant (A)',...
'$\quad  vf$ : volume fraction of the cylinder phase',...
'$\quad  drho$ : electron-density contrast',...
'$\quad  domainsize$ : crystal domain size (A)',...
'$\quad  microstrain$ : RMS strain broadening',...
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

Iq = hex_twophase(q, [p.d, p.vf], p.drho,p.domainsize,p.microstrain, p.DW);

back = p.poly1*q + p.poly2;
out = p.I0*Iq + back;

if isnan(out)
    out = ones(size(out))*1E20;
end

if nargout == 2
    a = p.d;
    volf = p.vf;
    
    area = a^2*sqrt(3)/2;
    R = sqrt(volf*area/pi);

    fprintf('Further information ======================================\n');
    fprintf('Lattice area : %0.3e %c^2.\n', area, char(197));
    fprintf('Radius of cylinder : %0.3e %c.\n', R, char(197));
    fprintf('==============================================================\n');
    
    report = '';
end
