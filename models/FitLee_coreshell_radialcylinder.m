function [out, report] = FitLee_coreshell_radialcylinder(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Polydisperse core-shell cylinder, radial direction fit in absolute unit $(cm^{-1}). $' ,...
'Data needs to be plotted as q vs $I(q) \cdot q$',...
'',...
'$I(q) = (volf/\left<V_{cyl}\right>) \cdot r_e^2 \cdot \left<|F_{radial}|^2\right> \cdot H^2 \cdot S(q) + I_b$',...
'$\qquad    S(q) = S_{HS}(q; D, v_f) + P_{cluster}(q; I_{0,2}, s_2, Rg_2, P_2)$',...
'$\qquad    I_b = poly1 \cdot q^{poly2} + poly3 \cdot q + poly4$',...
'',...
'$Parameters$',...
'$\quad  volf$ : volume fraction of cylinders (dimensionless)',...
'$\quad  coreR, shellthick$ : core radius and shell thickness (A)',...
'$\quad  coreR\_frac, shellthick\_frac$ : Schultz relative widths for core and shell',...
'$\quad  H$ : cylinder length (A); treated as fixed (no length polydispersity)',...
'$\quad  core\_edensity, shell\_edensity, solvent\_edensity$ : electron density in $A^{-3}$',...
'$\quad  D$ : interparticle distance (A) for hard-sphere $S(q)$',...
'$\quad  vf$ : volume fraction in hard-sphere $S(q)$',...
'$\quad  I0\_2, s\_2, Rg2, P2$ : Guinier-Porod cluster scatterer parameters',...
'',...
'Note: $r_e$ = 2.818E-5 A is the classical electron radius.',...
'$  I(q) units: (A^{-3}) * (A^2) * (A^0) = A^{-1}, /1E-8 (cm/A) -> cm^{-1}.$',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
};

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
    bestP.volf = 0.1;
    bestP.coreR = 25;
    bestP.coreR_frac = 0.1;
    bestP.shellthick = 10;
    bestP.shellthick_frac = 0.1;
    bestP.core_edensity = 0;            % electron density in A^-3
    bestP.shell_edensity = 0.39;
    bestP.solvent_edensity = 0.3344;
    bestP.H = 200;

    bestP.I0_2 = 1;
    bestP.s_2 = 0;
    bestP.Rg_2 = 300;
    bestP.P_2 = 4;
    
    bestP.D = 100;
    bestP.vf = 0.01;
    bestP.poly1 = 0;
    bestP.poly2 = -2;
    bestP.poly3 = 0;
    bestP.poly4 = 0;
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

% radial direction
q = q(:);

if p.coreR_frac > 0
    numpnt = 21;
    [nr, x] = schultzdist99(p.coreR, p.coreR*p.coreR_frac, numpnt);
else
    numpnt = 1;
    nr = 1;
    x = p.coreR;
end
if p.shellthick_frac > 0
    numpnt2 = 21;
    [nr2, x2] = schultzdist99(p.shellthick, p.shellthick*p.shellthick_frac, numpnt2);
else
    numpnt2 = 1;
    nr2 = 1;
    x2 = p.shellthick;
end

Iq = zeros(size(q));
S = 0;
for j=1:numpnt2
    for i=1:numpnt
        sh = x2(j);
        cR = x(i);
        R = cR + sh;
        radius = x(i)/p.coreR*[cR, cR + sh];
        edensity = [p.core_edensity, p.shell_edensity, p.solvent_edensity];
        Iq0 = saxscylinder_CS(q, radius, edensity, 0);
        S = S + pi*R^2*nr(i)*nr2(j);   % was (pi*R)^2: off-by-pi
        Iq = Iq + nr2(j)*nr(i)*Iq0*p.H^2;
    end
end
Vcyl = S * p.H;

r_e = 2.818E-5; % Angstrom
Angstrom2Centimeter = 1E-8;

Sq0 = strfactor2(q, p.D, p.vf);
Pcluster = guinierporodmodel(q, p.I0_2, p.P_2, p.Rg_2, p.s_2);
Sq = Sq0 + Pcluster;
%Sq = 1;

pnumberfraction = p.volf/(Vcyl);
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;

out = pnumberfraction*r_e^2*Iq.*Sq/Angstrom2Centimeter + back;
if isnan(out)
    out = ones(size(out))*1E20;
end

if nargout == 2
    maxR = p.coreR+p.shellthick;
    sigma_rel = max(p.coreR_frac, p.shellthick_frac);   % widest distribution for plot range
    x = linspace(maxR/10, maxR + maxR*sigma_rel*9, 100);
    nr = schultzdist(x, maxR, maxR*sigma_rel);

    figure;
    plot(x, nr);xlabel('Radius (A)');ylabel('n(R)')
    report = '';
end


