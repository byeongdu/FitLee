function [out, report] = FitLee_coreshellshell_polymer(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Polydisperse core-shell-shell sphere with polymer chains, absolute unit $(cm^{-1}). $' ,...
'',...
'$I(q) = r_e^2 \cdot \left\{ (fn0/\left<V_p\right>) \cdot \left( \left<|F_{cs}|^2\right> \cdot S_0(q) + N_c \cdot |F_{chain}|^2 \cdot S_1(q) \right) \right.$',...
'$\qquad\qquad \left. + (fn1/\left<V_1\right>) \cdot \Delta\rho_1^2 \cdot \left<|F_{sph}|^2\right> \cdot S_1(q) \right\} + I_b$',...
'$\qquad    S\_0(q) = S_{HS}(q; D0, vf0), \quad S\_1(q) = S_{HS}(q; D, vf)$',...
'$\qquad    I_b = poly1 \cdot q^{poly2} + poly3 \cdot q + poly4$',...
'',...
'$Parameters$',...
'$\quad  fn0$ : volume fraction of core-shell-shell particles (dimensionless)',...
'$\quad  coreRdivshR, sigCore$ : core-to-shell radius ratio and its rel. width',...
'$\quad  shellR, sigshellR$ : shell radius peak (A) and rel. width',...
'$\quad  sh2thick, sigsh2thick$ : 2nd shell thickness peak (A) and rel. width',...
'$\quad  core\_edensity, sh1\_edensity, sh2\_edensity, solvent\_edensity$ : electron density in $A^{-3}$',...
'$\quad  d\_edensity\_chain$ : chain-vs-solvent electron density contrast $(A^{-3})$',...
'$\quad  Nc, Rg$ : number of chains per particle and chain radius of gyration (A)',...
'$\quad  fn1, d\_edensity1, r1, sig1$ : optional 2nd sphere population $(fn1 = 0$ disables$)$',...
'$\quad  D0, vf0$ : $S(q)$ for core-shell-shell particles',...
'$\quad  D, vf$ : $S(q)$ for chains and 2nd sphere',...
'',...
'Note: $r_e$ = 2.818E-5 A is the classical electron radius.',...
'$  I(q) units: (A^{-3}) * (A^2) * (A^0) = A^{-1}, /1E-8 (cm/A) -> cm^{-1}.$',...
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
    bestP.D0 = 100;
    bestP.vf0 = 0.01;
    bestP.fn0 = 0.01;                 % volume fraction of core-shell-shell particles
    bestP.coreRdivshR = 0.245672;
    bestP.sigCore = 0.097654;
    bestP.shellR = 247.691904;
    bestP.sigshellR = 0.012843;
    bestP.sh2thick = 42.636334;
    bestP.sigsh2thick = 0.046388;
    bestP.core_edensity = 0.38;            % electron density in A^-3
    bestP.sh1_edensity = 0.39;
    bestP.sh2_edensity = 0.42;
    bestP.solvent_edensity = 0.3344;
    bestP.d_edensity_chain = 0.05;         % chain-vs-solvent electron density contrast (A^-3)
    bestP.Nc = 100;
    bestP.Rg = 50;
    % Optional second sphere population (disable with fn1 = 0)
    bestP.fn1 = 0;                         % volume fraction
    bestP.d_edensity1 = 0.05;              % electron density contrast (A^-3)
    bestP.r1 = 35;
    bestP.sig1 = 3.5;
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

method = 1;

% Polymer chain contribution (per chain): (d_edensity*V_chain)^2 * Debye(qRg)  [units A^0]
Rg = p.Rg;
Vchain = 4*pi/3*Rg^3;
Pc_perchain = (p.d_edensity_chain*Vchain)^2 * Pchain(q, Rg);
Pc_perparticle = p.Nc * Pc_perchain;     % Nc chains per core-shell-shell particle

% Core-shell-shell: <|F_cs|^2> and matching <V_p> from the SAME distribution grid
[Iq_cs, V_p_mean] = coreshell(p, q, method);

% Structure factors and background
Sq0 = strfactor2(q, p.D0, p.vf0);
Sq1 = strfactor2(q, p.D, p.vf);
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;

% Absolute scale: I[cm^-1] = (fn0/<V_p>)[A^-3] * r_e^2[A^2] * <|F|^2>[A^0] / 1E-8 [cm/A]
pscale_cs = p.fn0 / V_p_mean * r_e^2 / Angstrom2Centimeter;
Iq1 = pscale_cs * Iq_cs(:)        .* Sq0;
Pc  = pscale_cs * Pc_perparticle(:) .* Sq1;

% Optional second sphere population
if p.fn1 > 0
    [Pq1_norm, V1_2, V1_1] = SchultzsphereFun(q, p.r1, p.sig1);
    Pq1_F2  = V1_2 * Pq1_norm(:);   % <|F|^2 / d_edensity^2> in A^6
    pscale1 = (p.fn1 / V1_1) * p.d_edensity1^2 * r_e^2 / Angstrom2Centimeter;
    Iq2 = pscale1 * Pq1_F2 .* Sq1;
else
    Iq2 = zeros(numel(q), 1);
end

out = Iq1 + Iq2 + Pc + back;
out = [out(:), Iq1(:), Iq2(:), Pc(:), back(:)];

if any(isnan(out(:)))
    out = ones(size(out));
end

if nargout == 2
    x = 0:1:(p.shellR + p.sh2thick + 10);
    rho = zeros(size(x));
    t = x < p.coreRdivshR*p.shellR;
    
    rho(t) = p.core_edensity;
    t = (x >= p.coreRdivshR*p.shellR) & (x < p.shellR);
    rho(t) = p.sh1_edensity;
    t = (x >= p.shellR) & (x < p.shellR+p.sh2thick);
    rho(t) = p.sh2_edensity;
    t = (x >= p.shellR+p.sh2thick);
    rho(t) = p.solvent_edensity;
    
    figure;
    plot(x, rho);xlabel('Radius (A)');ylabel('\rho (R)')
    report = '';
end

function y = Pchain(q, Rg)
    x = (q*Rg).^2;
    y = 2*(exp(-x)-1+x)./x.^2;

function [Iq, V_p_mean] = coreshell(p, q, method)
% V_p_mean : mean total particle volume <V_p> = <(4pi/3)*(shellR + sh2thick)^3>
%            computed from the SAME nr2, nr3 distribution grids used for <|F|^2>,
%            so the per-particle factor cancels in (fn0/<V_p>) * <|F|^2>.

V_p_mean = NaN;   % set inside method branches

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
    % <V_p> over the (r2, r3) grid; uses identical nr2/nr3 sums to the loop above,
    % so any truncation/normalization drift cancels in (fn0/<V_p>) * <|F|^2>.
    Vp_acc = 0; norm_acc = 0;
    for mm = 1:numel(nr2)
        for ll = 1:numel(nr3)
            w = nr2(mm) * nr3(ll);
            R_outer = r2(mm) + r3(ll);
            Vp_acc   = Vp_acc   + w * (4*pi/3) * R_outer^3;
            norm_acc = norm_acc + w;
        end
    end
    V_p_mean = Vp_acc / norm_acc;
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