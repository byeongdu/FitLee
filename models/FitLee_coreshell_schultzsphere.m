function [out, report] = FitLee_coreshell_schultzsphere(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Polydisperse core-shell sphere fit in absolute unit $(cm^{-1}). $' ,...
'',...
'$I(q) = (fn0/\left<V_p\right>) \cdot r_e^2 \cdot \left<|F(q)|^2\right> \cdot S(q; D, v_f) + I_b$',...
'$\qquad    \left<|F(q)|^2\right>$ from SchultzCoreShellFF3 with electron-density-weighted amplitude',...
'$\qquad    I_b = poly1\cdot q^{poly2} + poly3 \cdot q + poly4$',...
'',...
'$Parameters$',...
'$\quad  fn0$ : volume fraction of particles (dimensionless)',...
'$\quad  coreR$ : core radius peak (A)',...
'$\quad  sigma$ : size distribution width (lognormal sigma, dimensionless)',...
'$\quad  shellthick$ : shell thickness (A); $Rshell = coreR + shellthick$ scales with $R$',...
'$\quad  core\_edensity, shell\_edensity, solvent\_edensity$ : electron density in $A^{-3}$',...
'$\quad  D$ : interparticle distance (A) for hard-sphere $S(q)$',...
'$\quad  vf$ : volume fraction in hard-sphere $S(q)$',...
'',...
'Note: $r_e$ = 2.818E-5 A is the classical electron radius.',...
'$  I(q) units: (A^{-3}) * (A{^2}) * (A^0) = A^{-1}, /1E-8 (cm/A) -> cm^{-1}.$',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    1. Kwon et al. Nature Materials, 2015, 14, 215-223. ',...
'    2. Wang et al. J. Phys. Chem. C. 2013. 117(44), 22627. '};

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
    bestP.fn0 = 0.01; % volume fraction of particles (dimensionless)
    bestP.coreR = 25;
    bestP.sigma = 0.1;
    bestP.shellthick = 10;
    bestP.core_edensity = 0;
    bestP.shell_edensity = 0.39;
    bestP.solvent_edensity = 0.3344;
    bestP.D = 100;
    bestP.vf = 0.0;
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

radius = [p.coreR, p.coreR + p.shellthick];
edensity = [p.core_edensity, p.shell_edensity, p.solvent_edensity];

% <|F(q) * d_edensity|^2> averaged over polydisperse core radius
% Units: (A^-3)^2 * (A^3)^2 = A^0 (dimensionless); r_e^2 supplies the A^2.
Iq = SchultzCoreShellFF3(q, [radius, p.sigma, edensity]);

% Mean total particle volume <V_p> from the SAME distribution
% SchultzCoreShellFF3 uses internally (logndist9 on coreR, with
% Rshell = beta*Rcore where beta = (coreR+shellthick)/coreR).
[gDist, Rcore_vec] = logndist9(p.coreR, p.sigma);
beta = (p.coreR + p.shellthick) / p.coreR;
V_p_mean = sum(gDist .* (4*pi/3 .* (beta*Rcore_vec).^3));  % A^3

% Number density from volume fraction: N/V_sample = fn0/<V_p>
pnumberfraction = p.fn0 / V_p_mean;   % A^-3

Sq1 = strfactor2(q, p.D, p.vf);
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;

% pnumberfraction[A^-3] * r_e^2[A^2] * Iq[A^0] = A^-1; /1E-8 (cm/A) -> cm^-1
out = pnumberfraction * r_e^2 * Iq(:) .* Sq1(:) / Angstrom2Centimeter + back;

if any(isnan(out))
    out = ones(size(out));
end

if nargout == 2
    maxR = p.coreR + p.shellthick;
    x = linspace(maxR/10, maxR + maxR*p.sigma*9, 200);
    nr = logndist(p.coreR, p.sigma, x);   % core radius distribution

    figure;
    plot(x, nr); xlabel('Core radius (A)'); ylabel('n(R)')

    Rcore_mean    = sum(gDist .* Rcore_vec);
    Rshell_outer  = beta * Rcore_mean;
    pnumdens_cm3  = pnumberfraction / Angstrom2Centimeter^3;  % #/cm^3
    fprintf('Statistical information of the core-shell particle ===================\n');
    fprintf('Mean core radius          : %0.3e %c.\n',   Rcore_mean, char(197));
    fprintf('Mean outer radius         : %0.3e %c.\n',   Rshell_outer, char(197));
    fprintf('Mean particle volume <V_p>: %0.3e %c^3.\n', V_p_mean, char(197));
    fprintf('Number concentration      : %0.3e particles/cm^3.\n', pnumdens_cm3);
    fprintf('Mol concentration         : %0.3e M (mole/L).\n',     pnumdens_cm3/6.022E23/1E-3);
    fprintf('Weight concentration (g/mL) = fn0 * particle density (g/mL).\n');
    fprintf('======================================================================\n');
    report = '';
end


function generate_Rmatrix(minR, maxR, numR, q)

    R = linspace(minR, maxR, numR);
    Rmat.r1 = R;
    Rmat.Fmat = sphereamp(q, R);
    Rmat.qmin = min(q);
    Rmat.qmax = max(q);
    assignin('base', 'Rmat', Rmat)

function Iq = coreshell(p, q, method)

% method = 1;

if method == 1
    if p.sigCore ==0
        r1 = p.coreRdivshR;
        nr1 = 1;
    else
        r1 = linspace(0.01, 0.99, 15);
        %nr1 = schultz(p.coreRdivshR, p.sigCore, r1);
        nr1 = schultz(p.coreRdivshR, p.sigCore*p.coreRdivshR, r1);
    end
    if p.sigshellR==0
        r2 = p.shellR;
        nr2 = 1;
    else
        %r2 = linspace(1, p.shellR+3*p.sigshellR, 30);
        r2 = linspace(1, p.shellR+3*p.sigshellR*p.shellR, 30);
        %nr2 = schultz(p.shellR, p.sigshellR, r2);
        nr2 = schultz(p.shellR, p.sigshellR*p.shellR, r2);
    %Rmat.r2 = r2;
    end
%     Rmat.r2 = r2;
%     Rmat.dist2 = nr2;
    if p.sigsh2thick ==0
        r3 = p.sh2thick;
        nr3 = 1;
    else
%         r3 = linspace(0.1, p.sh2thick+3*p.sigsh2thick, 30);
%         nr3 = schultz(p.sh2thick, p.sigsh2thick, r3);
        r3 = linspace(0.1, p.sh2thick+3*p.sigsh2thick*p.sh2thick, 30);
        nr3 = schultz(p.sh2thick, p.sigsh2thick*p.sh2thick, r3);
    end
%     Rmat.r3 = r3;
%     Rmat.dist3 = nr3;

    yf = zeros(size(q));
    for mm=1:numel(nr2);
        for kk=1:numel(nr1)
            rad = [];
            for ll=1:numel(nr3)
                if p.coreRdivshR == 0
                    rad = [rad;r2(mm),r2(mm) + r3(ll)];
                else
                    rad = [rad;r1(kk)*r2(mm), r2(mm), r2(mm)+r3(ll)];
                end
            end
            if p.coreRdivshR == 0
                eden = [p.sh1_edensity, p.sh2_edensity, p.solvent_edensity];
            else
                eden = [p.core_edensity, p.sh1_edensity, p.sh2_edensity, p.solvent_edensity];
            end
            yt = multilayersphere2(q, rad,eden, nr3);
            yf = yf + yt(:)*nr2(mm)*nr1(kk);
        end
    end
end
if method == 2
    Rmat = evalin('base', 'Rmat');

    if p.coreRdivshR ==0
        nr1 = 1;
        nr2 = schultz(p.shellR, p.sigshellR, r2);
        %nr2 = 1;
        %Rmat.dist2 = 1;
        %Rmat.dist2 = nr2;
        if p.sigsh2thick ==0
            nr3 = 1;
        else
            nr3 = schultz(p.shellR+p.sh2thick, p.sigsh2thick, r2);
        end
        %dist3 = nr3;
        eden = [p.sh1_edensity, p.sh2_edensity, p.solvent_edensity];
        dis = {nr2, nr3};
        Fm = {Rmat.Fmat2, Rmat.Fmat3};
    else
        eden = [p.core_edensity, p.sh1_edensity, p.sh2_edensity, p.solvent_edensity];
        nr1 = schultz(p.coreRdivshR*p.shellR, p.coreRdivshR*p.sigCore*p.shellR, Rmat.r1);
        Rmat.dist1 = nr1;
        nr2 = schultz(p.shellR, p.sigshellR, Rmat.r2);
        Rmat.dist2 = nr2;
        nr3 = schultz(p.sh2thick+p.shellR, p.sigsh2thick, Rmat.r2);
        Rmat.dist3 = nr3;
        dis = {nr1, nr2, nr3};
        Fm = {Rmat.Fmat1, Rmat.Fmat2, Rmat.Fmat3};
    end
    yf = multisphere(Fm, eden, dis);
    assignin('base', 'Rmat', Rmat)
end

Iq = yf(:)/sum(nr1)/sum(nr2)/sum(nr3);

function Pq = multisphere(Fm, eden, dis)
eden = eden(1:end-1) - eden(2:end);
Fmat1 = Fm{1};dist1 = dis{1};
Fmat2 = Fm{2};dist2 = dis{2};
if numel(Fm) == 3
    Fmat3 = Fm{3};
    dist3 = dis{3};
end
    
numq = size(Fmat1, 1);
if numel(Fm) == 2
    for i=1:numq
        [f1,f2]=meshgrid(Fmat1(i,:),Fmat2(i,:));
        [d1,d2] = meshgrid(dist1, dist2);
        f = f1*eden(1)+f2*eden(2);

        f = abs(f).^2.*d1.*d2;
        f = sum(sum(f));
        Pq(i) = f;
    end
end
if numel(Fm) == 3
    for i=1:numq
        [f1,f2,f3]=meshgrid(Fmat1(i,:),Fmat2(i,:),Fmat3(i,:));
        [d1,d2,d3] = meshgrid(dist1, dist2, dist3);
        f = f1*eden(1)+f2*eden(2)+f3*eden(3);
        f = abs(f).^2.*d1.*d2.*d3;
        f = sum(sum(sum(f)));
        Pq(i) = f;
    end
end