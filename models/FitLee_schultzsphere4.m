function [out, report] = FitLee_schultzsphere4(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Schultz polydisperse sphere, two populations with separate $S(q)$, absolute unit $(cm^{-1}). $' ,...
'',...
'$I(q) = (fn0/\left<V_0\right>) \cdot r_e^2 \cdot \left( \delta\rho_0^2 \cdot P_0(q) \cdot S_a(q) + N_{ratio} \cdot \delta\rho_1^2 \cdot P_1(q) \cdot S_1(q) \right) + I_b$',...
'$\qquad    S_a(q) = powI \cdot q^{PorodExp} + S_{HS}(q; D_0, v_{f,0})$',...
'$\qquad    S_1(q) = S_{HS}(q; D_1, v_{f,1})$',...
'$\qquad    I_b = poly1 \cdot q^{poly2} + poly3 \cdot q + poly4$',...
'',...
'Note: the $powI \cdot q^{PorodExp}$ low-q tail is added only to particle 0''s $S(q)$.',...
'Particle 1 uses $\left<V_0\right>$ (particle 0''s volume) for the number-density denominator.',...
'Output is convolved through apply_resolution before return.',...
'If $v_f \geq 1$, strfactor3 is used instead of strfactor2.',...
'',...
'$Parameters$',...
'$\quad  fn0$ : volume fraction of particle 0 (dimensionless)',...
'$\quad  delta\_rho0, delta\_rho1$ : electron density contrast $(A^{-3})$',...
'$\quad  r0, sigfrac0$ : Schultz radius peak (A) and relative width for particle 0',...
'$\quad  r1, sigfrac1$ : same for particle 1 (Schultz FWHM = $r \cdot sigfrac$)',...
'$\quad  Nratio$ : number ratio of particle 1 over particle 0',...
'$\quad  powI, PorodExp$ : low-q power-law term added into $S\_a(q)$ only',...
'$\quad  D0, vf0$ : hard-sphere $S_{HS}$ parameters for particle 0',...
'$\quad  D1, vf1$ : hard-sphere $S_{HS}$ parameters for particle 1',...
'',...
'Note: $r_e$ = 2.818E-5 A is the classical electron radius.',...
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
    Nf = p;
    bestP = [];
    bestP.delta_rho0 = 1;
    bestP.r0 = 100;
    bestP.sigfrac0 = 10;
    bestP.fn0 = 1;
    bestP.D0 = 100;
    bestP.vf0 = 0.0;
    bestP.powI = 1E-4;
    bestP.PorodExp = -4;
    bestP.delta_rho1 = 1;
    bestP.r1 = 100;
    bestP.sigfrac1 = 10;
    bestP.Nratio = 0;
    bestP.D1 = 100;
    bestP.vf1 = 0.0;
    % Need 4 parameters for background.
    bestP.poly1 = 0;
    bestP.poly2 = -2;
    bestP.poly3 = 0;
    bestP.poly4 = 0;
    assignin('base', 'bestP', bestP);
    out = bestP;
    return
end

%% fitting codes .......................
if iscell(q);
    if numel(q) > 1
        error('FitLee_schultzsphere.m is for fitting a set of data, for now')
    end
    q = q{1};
end

r_e = 2.818E-5; % Angstrom
Angstrom2Centimeter = 1E-8;
%r_e2 = (r_e*Angstrom2Centimeter)^2;

q = q(:);

[Pq1, V1, V0] = SchultzsphereFun(q, p.r0, p.sigfrac0*p.r0);
Pq1 = V1*Pq1(:);
if abs(p.Nratio) > 0
    [Pq2, V2] = SchultzsphereFun(q, p.r1, p.sigfrac1*p.r1);
    Pq2 = V2*Pq2(:);
else
    Pq2 = 0;
end

if p.vf0<1
    Sq0 = strfactor2(q, p.D0, p.vf0);
else
    Sq0 = strfactor3(q, [p.D0, p.vf0]);
end
if p.vf1<1
    Sq1 = strfactor2(q, p.D1, p.vf1);
else
    Sq1 = strfactor3(q, [p.D1, p.vf1]);
end
Sq = p.powI*q.^p.PorodExp + Sq0;
%Sq = p(4)*q(:).^p(5) + strfactor_2Dpara(q, p(6), p(7));
%Iq = p(1)*Imat*nr1/sum(nr1)
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
pnumberfraction = p.fn0/(V0);
Iq1 = pnumberfraction/Angstrom2Centimeter*r_e^2*p.delta_rho0^2*Pq1.*Sq;
Iq2 = pnumberfraction/Angstrom2Centimeter*r_e^2*p.Nratio*p.delta_rho1^2*Pq2.*Sq1;
%out = pnumberfraction*r_e^2*(p.delta_rho0^2*Pq1.*Sq+ p.Nratio*p.delta_rho1^2*Pq2.*Sq1);
out = Iq1 + Iq2 +back; % data is in cm^-1
out = apply_resolution(q, out);
out = [out(:), Iq1(:), Iq2(:), back(:)];
if isnan(out)
    out = ones(size(out));
end

if nargout == 2
    x = 0:1:(max(p.r0, p.r1)+max(p.sig0, p.sig1)*10);
    nr0 = schultz(p.r0, p.sig0, x);
    nr = nr0;
    if p.Nratio > 0
        nr1 = schultz(p.r1, p.sig1, x);
        nr = nr+p.Nratio*nr1;
    end
    figure;
    plot(x, nr);xlabel('Radius (A)');ylabel('n(r)')
    fprintf('Volume of particle0 : %0.3e A^3\n', V0);
    pnumberfraction = pnumberfraction/Angstrom2Centimeter^3; % num fraction in 1/cm^3 unit.

    fprintf('Number fraction of particle0 : %0.3e particles/cm^3. \n', pnumberfraction)
    fprintf('Mol concentration : %0.3e M (mole/L). \n', pnumberfraction/6.022E23/1E-3)
    ns = input('Type the  density of a particle in g/mL unit. e.g. water = 1;');
    fprintf('Weight concentration of your particles: %0.3e g/mL. \n', p.fn0*ns)
    
    report = '';
end