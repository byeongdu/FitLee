function [out, report] = FitLee_poly_fractalparticle2(varargin)
FitLee_helpstr = {'Polydisperse mass-fractal particles. ' ,...
'I(q) = I0 \cdot P(q; avgR, z, D_f) \cdot S(q) + I_b(q)',...
'    S(q) = guinierporodmodel(q, I0_2, s_2, Rg_2, P_2) + S(q; Rh, \nu)',...
'    I_b(q) = poly1 \cdot q^{poly2} + poly3 \cdot q + poly4',...
'Parameters',...
'\bfP(q): Form factor of PP',...
'  \rmI0 : Scale factor I0',...
'  avgR : Schultz size distribution, radius peak (A)',...
'  z: Schultz size distribution, FWHM (A)',...
'  D_f : mass fractal dimension. Set D_f<3 or 4.',...
'',...
'\bfS(q) = P_{cluster}(q) + S_{local}(q):',...
' \bfP_{cluster}(q) = guinierporodmodel(q, I0, d, Rg, s)',...
'   \rmI0 : Guinier I0, meaning the number of primary particles forming a cluster',...
'   Rg : Rg of the cluster',...
'   P : Porod exponent of the cluster',...
'   s : dimensionality of the cluster',... 
' \bfS_{local}(q)',...
'  \rmRh : Hydrodynamic radius (A) (hard sphere potential S(q))',...
'  \nu : volume fraction of particle 0 (hard sphere potential S(q)) ',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    1. R. Besselink and J. E. ten Elshof, J. Appl. Cryst., (2014), 47, 1606. '};

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
    bestP.I0 = 1;
    bestP.r1 = 100;
    bestP.sig1 = 10;
    bestP.Rh1 = 100;
    bestP.vf1 = 0.0;
    % structure factor
    bestP.Nratio2 = 1;
    bestP.s2 = 0;
    bestP.Rg2 = 300;
    bestP.P2 = 4;
    bestP.Rh2 = 100;
    bestP.vf2 = 0.0;
%    bestP.powI_Cluster = 1E-8;
%    bestP.PorodExp_Cluster = -4;
    % Need 4 parameters for background.
    bestP.poly1 = 0;
    bestP.poly2 = -2;
    bestP.poly3 = 0;
    bestP.poly4 = 0;
    bestP.poly5 = 0;
    assignin('base', 'bestP', bestP);
    out = bestP;
    return
end

%% fitting codes .......................
r_e = 2.818E-5; % Angstrom
Angstrom2Centimeter = 1E-8;

if iscell(q)
    if numel(q) > 1
        error('FitLee_schultzsphere.m is for fitting a set of data, for now')
    end
    q = q{1};
end

q = q(:);

try
    UserBack = evalin('base', 'FitLeeUserBack');
    UserBack = interp1(UserBack(:,1),UserBack(:,2),q);
catch
    UserBack = 0;
end

RgPP = [];
[Pq1, V2] = SchultzsphereFun(q, p.r1, p.sig1);
Pq1 = V2*Pq1;

pnumberfraction = p.I0/(sqrt(V2));
p.I0 = pnumberfraction*r_e^2/Angstrom2Centimeter;

Sq2 = strfactor2(q, p.Rh2, p.vf2);
Sq1 = strfactor2(q, p.Rh1, p.vf1);
Pq2 = guinierporodmodel(q, p.Nratio2, p.P2, p.Rg2, p.s2);
Pq2 = Pq2*p.Rg2^6;
%Poq = p.powI_Cluster*q.^p.PorodExp_Cluster;
Iq1 = Pq1.*Sq1;
Iq2 = p.Nratio2*Pq2.*Sq2;
Iq = p.I0*(Iq1 + Iq2);
back = p.poly1*q(:).^p.poly2 + p.poly3*q(:) + p.poly4 + p.poly5*UserBack;
%pnumberfraction = p.f0;

out = [Iq + back, p.I0*Iq1(:), p.I0*Iq2(:), back];
if isnan(out)
    out = ones(size(out));
end

if nargout == 2
    x = 0:1:((p.r0_PP+p.sig0_PP)*10);
    nr0 = schultz(p.r0_PP, p.sig0_PP, x);
    nr = nr0;
    figure;subplot(2,1,1)
    plot(x, nr);xlabel('Radius (A)');ylabel('n(r)')
    subplot(2,1,2)
    Vr = nr(:).*x(:).^3;
    plot(x, Vr);xlabel('Radius (A)');ylabel('V(r)')
    volumedistribution = [x(:), Vr];
    assignin('base', 'volumedistribution', volumedistribution);
    %[zRg, V, ~, S] = schultzRg(p.r0, p.sig0);
    if ~isempty(RgPP)
        zRg = RgPP;
    else
        zRg = 0;
    end
    fprintf('Rg of particle0 : %0.3fA.\n', zRg);
    %[~, ~, ~, S] = schultzRg(p.r0, p.sig0);
    %S = S/4/100; % 4 pi R(A) ^2 --> pi*R(nm)^2
    report = '';
end
