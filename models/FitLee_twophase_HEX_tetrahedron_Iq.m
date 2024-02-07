function [out, report] = FitLee_twophase_HEX_tetrahedron_Iq(varargin)
    % Paracrystal theory equation for triblock copolymer.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'In-plane diffraction of 2D hexagonally packing of tetrahedrons.' ,...
'Tetrahedrons are at the lattice point of 2D hexagon.',...
'In a unit cell, there is one tetrahedron at the lattice points.',...
'Its base triangle is parallel to axb plane.',...
'',...
'Note that this function calculates up to [5,0,0] reflections.',...
'Therefore, qmax should be 2*pi/(a*sqrt(3)/2)*5.',...
'',...
'2D hexagonally packed tetrahedrons',...
'I0 : Scalling constant',...
'a : Lattice Constant (in Angstrom)',...
'L : Edge length of a tetrahedron (A)',...
'th0 : Orientation angle of the tetrahedron (degree)',...
'gaussian : Gaussian width of the Voight diffraction peak',...
'domainsize : Crystal domain size that defines the Lorentz width',...
'DW : Debye-Waller factor',...
'',...
'    background = poly1*q + poly2',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
' Senesi, A. J. and B. Lee (2015). ',...
' Small-Angle Scattering of Particle Assemblies. JAC. 48(4): 1172-1182.'};

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
    bestP.a = 200; % Angstrom
    bestP.L = 190; % Angstrom
    bestP.th0 = 0; % degree
    bestP.domainsize = 5000;
    bestP.gaussian = 0.0025;
    bestP.microstrain = 0.02;
    bestP.DW = 0.02;
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

Sq = hex_tetrahedron(q, [p.a, p.L, p.th0], p.domainsize,p.gaussian, p.microstrain, p.DW);
Pq = hex_tetrahedron(q, [p.a, p.L, p.th0], 'pq', true);

back = p.poly1*q + p.poly2;
out = p.I0*Sq.*Pq + back;

if isnan(out)
    out = ones(size(out))*1E20;
end

if nargout == 2
    a = p.a;
    L = p.L;

    fprintf('Further information ======================================\n');
    fprintf('Lattice constant : %0.3e %c.\n', a, char(197));
    fprintf('Edge length : %0.3e %c.\n', L, char(197));
    fprintf('==============================================================\n');

    report = '';
end
