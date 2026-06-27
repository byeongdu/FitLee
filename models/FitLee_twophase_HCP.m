function [out, report] = FitLee_twophase_HCP(varargin)
    % Paracrystal theory equation for triblock copolymer.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Two-phase 3D HCP lattice diffraction.' ,...
'',...
'$I(q) = I_0 \cdot $ hcp_twophase$(q, [d, v_f], \Delta\rho, domainsize, microstrain, DW) + I_b$',...
'$\qquad    I_b = poly1 \cdot q + poly2$',...
'',...
'HCP unit cell: 2 spheres of phase A; lattice with $a = d$ and ideal $c = \sqrt{8/3} \cdot a$.',...
'Unit-cell volume $V = a^2 \cdot \sqrt{3}/2 \cdot c$; sphere radius $R = (v_f V / 2 / (4\pi/3))^{1/3}$.',...
'',...
'$Parameters$',...
'$\quad  I0$ : overall scale factor',...
'$\quad  d$ : lattice constant $a$ (A)',...
'$\quad  vf$ : volume fraction of phase A',...
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
Iq = hcp_twophase(q, [p.d, p.vf], p.drho,p.domainsize,p.microstrain, p.DW);

back = p.poly1*q + p.poly2;
out = p.I0*Iq + back;

if isnan(out)
    out = ones(size(out))*1E20;
end

if nargout == 2
    a = p.d;
    volf = p.vf;
    
    c = sqrt(8/3)*a;
    vol = a^2*sqrt(3)/2*c;
    
    R = (volf*vol/2/(4*pi/3))^(1/3);
    fprintf('Further information ======================================\n');
    fprintf('Lattice parameters : a = %0.3e and c = %0.3e %c.\n', a, c, char(197));
    fprintf('Lattice volume : %0.3e %c^3.\n', vol, char(197));
    fprintf('Radius of sphere : %0.3e %c.\n', R, char(197));
    fprintf('==============================================================\n');
    
    report = '';
end
