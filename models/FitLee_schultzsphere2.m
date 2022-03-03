function [out, report] = FitLee_schultzsphere2(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Schultz polydisperse sphere fit in absolute unit. ' ,...
'$I(q) = fn_0\cdot(Sq\cdotP(q; r_0, \sigma_0) + N_{ratio}\cdot\delta_{\rho1}^2*P(q; r_1, \sigma_1)) + Ib$',...
'    Sq = powI*q.^PorodExp + S(q; D, vf)',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'Parameters',...
'  fn0 : volume fraction of particle 0 (A^-3)',...
'  delta_rho0 : electron density difference for particle 0 (A^-3)',...
'  r0 : Schultz size distribution, radius peak (A)',...
'  sig0: Schultz size distribution, FWHM (A)',...
'  Nratio : Number ratio of particle 1 over particle 0',...
'  delta_rho1 : electron density difference for particle 1 (A)',...
'  D : Interparticle distance (A) (hard sphere potential S(q))',...
'  vf : volume fraction of particle 0 (hard sphere potential S(q)) ',...
'',...
'Note that',...
'I(q)_diff_cross_section = (delta_rho*r_e)^2*f_n*P(q)',...
'         , where P(0) = V_p^2'...
' ',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    1. Kwon et al. Nature Materials, 2015, 14, 215–223. ',...
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
    bestP.sig0 = 10;
    bestP.fn0 = 1;
    bestP.delta_rho1 = 1;
    bestP.r1 = 100;
    bestP.sig1 = 10;
    bestP.Nratio = 0;
    bestP.D = 100;
    bestP.vf = 0.0;
    bestP.powI = 1E-4;
    bestP.PorodExp = -4;
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
if iscell(q)
    if numel(q) > 1
        error('FitLee_schultzsphere.m is for fitting a set of data, for now')
    end
    q = q{1};
end

r_e = 2.818E-5; % Angstrom
Angstrom2Centimeter = 1E-8;
%r_e2 = (r_e*Angstrom2Centimeter)^2;

q = q(:);

[Pq1, V1, V0] = SchultzsphereFun(q, p.r0, p.sig0);
Pq1 = V1*Pq1(:);
if abs(p.Nratio) > 0
    [Pq2, V2] = SchultzsphereFun(q, p.r1, p.sig1);
    Pq2 = V2*Pq2(:);
else
    Pq2 = 0;
end

Sq1 = strfactor2(q, p.D, p.vf);
Sq = p.powI*q.^p.PorodExp + Sq1;
%Sq = p(4)*q(:).^p(5) + strfactor_2Dpara(q, p(6), p(7));
%Iq = p(1)*Imat*nr1/sum(nr1)
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
pnumberfraction = p.fn0/(V0*Angstrom2Centimeter^3);
out = pnumberfraction*r_e^2*Angstrom2Centimeter^2*(p.delta_rho0^2*Pq1.*Sq+ p.Nratio*p.delta_rho1^2*Pq2)+back;
if isnan(out)
    out = ones(size(out));
end

if nargout == 2
    x = 0:1:(max(p.r0, p.r1)+max(p.sig0, p.sig1)*10);
    x = x(:);
    nr0 = schultz(p.r0, p.sig0, x);
    nr = nr0;
    if p.Nratio > 0
        nr1 = schultz(p.r1, p.sig1, x);
        nr = nr+p.Nratio*nr1;
    end
    [zRg, V, V2, S] = schultzRg(p.r0, p.sig0);
    vmeanR = trapz(x, x.^3.*nr0);vmeanR = vmeanR^(1/3);
    v2meanR = trapz(x, x.^6.*nr0);v2meanR = v2meanR^(1/6);
    figure;
    plot(x, nr);xlabel('Radius (A)');ylabel('n(r)')
    fprintf('Statistical information of the particle 0 ======================================\n');
    fprintf('Number-mean radius of a single particle : %0.3e %c.\n', p.r0, char(197));
    fprintf('Volume-mean radius of a single particle : %0.3e %c.\n', vmeanR, char(197));
    fprintf('Volume-squared-mean radius of a single particle : %0.3e %c.\n', v2meanR, char(197));
    fprintf('Mean surface area of a single particle : %0.3e %c^2.\n', S, char(197));
    fprintf('Mean volume of a single particle : %0.3e %c^3.\n', V, char(197));
    fprintf('Mean square volume of a single particle : %0.3e %c^6.\n', V2, char(197));
    fprintf('Radius of gyration (Rg) : %0.3e %c.\n', zRg, char(197));
    fprintf('Number concentration : %0.3e particles/cm^3. \n', pnumberfraction)
    fprintf('Mol concentration : %0.3e M (mole/L). \n', pnumberfraction/6.022E23/1E-3)
    fprintf('Weight concentration (g/mL) can be obtained by multiplying your particles'' density (g/mL) to the fn0.\n');
    fprintf('==============================================================\n');
%     ns = input('Type the  density of a particle in g/mL unit. e.g. water = 1;');
%     fprintf('Weight concentration of your particles: %0.3e g/mL. \n', p.fn0*ns)
    
    report = '';
end
