function [out, report] = FitLee_coreshell_cylinder(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Polydisperse core-shell cylinder fit in absolute unit $(cm^{-1}). $' ,...
'Data needs to be plotted as q vs $I(q) \cdot q$ (radial integration)',...
'',...
'$I(q) = (volf/\left<V_{cyl}\right>) \cdot r_e^2 \cdot \left<|F_{radial}|^2 \cdot |F_{length}|^2\right> \cdot S(q; D, v_f) + I_b$',...
'$\qquad    I_b = poly1 \cdot q^{poly2} + poly3 \cdot q + poly4$',...
'',...
'$Parameters$',...
'$\quad  volf$ : volume fraction of cylinders (dimensionless)',...
'$\quad  coreR, shellthick$ : core radius and shell thickness (A); both scale with size',...
'$\quad  sigma$ : radial size distribution width (Schultz, relative)',...
'$\quad  H, sigH$ : cylinder length peak (A) and relative width (Schultz)',...
'$\quad  core\_edensity, shell\_edensity, solvent\_edensity$ : electron density in $A^{-3}$',...
'$\quad  D$ : interparticle distance (A) for hard-sphere $S(q)$',...
'$\quad  vf$ : volume fraction in hard-sphere $S(q)$',...
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
    bestP.volf = 0.01;                 % volume fraction of cylinders (dimensionless)
    bestP.coreR = 25;
    bestP.sigma = 0.1;
    bestP.shellthick = 10;
    bestP.core_edensity = 0;
    bestP.shell_edensity = 0.39;
    bestP.solvent_edensity = 0.3344;
    bestP.H = 200;
    bestP.sigH = 0.1;

    bestP.I0_2 = 1;
    bestP.s_2 = 2;
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

q = q(:);
Angstrom2Centimeter = 1E-8;
r_e = 2.818E-5;   % classical electron radius (A)
eden = [p.core_edensity, p.shell_edensity, p.solvent_edensity];   % A^-3

% --- Radial distribution (outer radius polydispersity; core/shell scale together) ---
maxR = p.coreR + p.shellthick;
if p.sigma > 0
    numpntR = 41;
    xR = linspace(sqrt(maxR/10), sqrt(maxR + maxR*p.sigma*9), numpntR).^2;
    nrR = schultzdist(xR, maxR, maxR*p.sigma);
else
    xR = maxR;  nrR = 1;
end
xR = xR(:);  nrR = nrR(:);
beta_core = p.coreR / maxR;   % preserve core/outer ratio across the distribution

% --- Length distribution ---
if p.H == 0
    xH = 0;  nrH = 1;
elseif p.sigH > 0
    numpntH = 55;
    xH = linspace((p.H/10).^(1/1.5), (p.H + p.H*p.sigH*9).^(1/1.5), numpntH).^1.5;
    nrH = schultzdist(xH, p.H, p.H*p.sigH);
else
    xH = p.H;  nrH = 1;
end
xH = xH(:).';  nrH = nrH(:).';

% --- Double loop: accumulate <|F|^2> and <V_cyl> with identical distribution weights ---
% Same sum(nrR)*sum(nrH) factor appears in both, so volf/V_acc * Iq cancels it exactly.
Iq = zeros(numel(q), 1);
V_acc = 0;
for j = 1:numel(nrH)
    H = xH(j);
    if H == 0
        Fl2 = ones(numel(q), 1);   % treat as cross-section only
        Vlen = 1;
    else
        Fl2 = sinc(q*H/2).^2 * H^2;   % preserved from original (note: MATLAB sinc convention)
        Vlen = H;
    end
    for i = 1:numel(nrR)
        R_outer = xR(i);
        radii   = [beta_core*R_outer, R_outer];
        Fr2 = saxscylinder_CS(q, radii, eden, 0);  % |F_radial|^2 (single radius)
        w   = nrR(i) * nrH(j);
        Iq    = Iq    + w * Fr2 .* Fl2;
        V_acc = V_acc + w * pi * R_outer^2 * Vlen;
    end
end

% Absolute scale: I(q) = (volf/<V_cyl>) * r_e^2 * <|F|^2> ; both <> share the same nr-sums.
% Units: A^-3 * A^2 * A^0 = A^-1, then /1E-8 (cm/A) -> cm^-1.
pscale = p.volf / V_acc * r_e^2 / Angstrom2Centimeter;

Pq_1 = guinierporodmodel(q, p.I0_2, p.P_2, p.Rg_2, p.s_2);
Sq1  = strfactor2(q, p.D, p.vf) + Pq_1;
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;

out = pscale * Iq(:) .* Sq1(:) + back;
if any(isnan(out))
    out = ones(size(out));
end

if nargout == 2
    maxR = p.coreR + p.shellthick;
    xplot = linspace(maxR/10, maxR + maxR*p.sigma*9, 100);
    nrplot = schultzdist(xplot, maxR, maxR*p.sigma);

    V_cyl_mean = V_acc / (sum(nrR)*sum(nrH));   % undo the common nr-sum for display
    pnumdens_cm3 = (p.volf / V_cyl_mean) / Angstrom2Centimeter^3;
    figure;
    plot(xplot, nrplot); xlabel('Outer radius (A)'); ylabel('n(R)')
    fprintf('Core-shell cylinder ==================================================\n');
    fprintf('Mean outer radius           : %0.3e %c.\n',  maxR, char(197));
    fprintf('Mean length                 : %0.3e %c.\n',  p.H,  char(197));
    fprintf('Mean cylinder volume <V_cyl>: %0.3e %c^3.\n', V_cyl_mean, char(197));
    fprintf('Number concentration        : %0.3e particles/cm^3.\n', pnumdens_cm3);
    fprintf('Mol concentration           : %0.3e M.\n',   pnumdens_cm3/6.022E23/1E-3);
    fprintf('======================================================================\n');
    report = '';
end


