function [out, report] = FitLee_EllipsoidalCoreShell(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Schultz polydisperse sphere fit. ' ,...
'I(q) = I0*(Sq*P(q; r0, sig0) + I1*P(q; r1, sig1)) + background',...
'    Sq = powI*q.^PorodExp + S(q; D, vf)',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'CF. if you have data in absolute unit, consider use FitLee_schultzsphere2',...
'',...
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
    bestP.I0 = 4.6095e-13;
    bestP.a = 30;
    bestP.eccentricity = 0.25;
    bestP.shell = 10;
    bestP.sig_a = 0.1;
    bestP.sig_shell = 0.0;
    bestP.core_rho = 0.39;
    bestP.shell_rho = 0.26;
    bestP.solvent_rho = 0.3344;
    bestP.I1 = 0;
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
if iscell(q);
    if numel(q) > 1
        error('FitLee_schultzsphere.m is for fitting a set of data, for now')
    end
    q = q{1};
end

q = q(:);

method = 1;

Iq = coreshell(p, q, method);

[Pq1, V1] = SchultzsphereFun(q, p.r1, p.sig1);
Pq1 = V1*Pq1(:);
Sq1 = strfactor2(q, p.D, p.vf);
%Sq = p.powI*q.^p.PorodExp + Sq1;

back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
I_ellip = p.I0*Iq.*Sq1;
I_sph = p.I1*Pq1;
Iq = I_ellip + I_sph + back;
out = [Iq(:), I_ellip(:), I_sph(:), back(:)];
if isnan(out)
    out = ones(size(out))*1E20;
end

if nargout == 2
    x = 0:1:(p.shellR + p.sh2thick + 10);
    rho = zeros(size(x));
    t = x < p.coreRdivshR*p.shellR;
    
    rho(t) = p.core_rho;
    t = (x >= p.coreRdivshR*p.shellR) & (x < p.shellR);
    rho(t) = p.sh1_rho;
    t = (x >= p.shellR) & (x < p.shellR+p.sh2thick);
    rho(t) = p.sh2_rho;
    t = (x >= p.shellR+p.sh2thick);
    rho(t) = p.solvent_rho;
    
    figure;
    plot(x, rho);xlabel('Radius (A)');ylabel('\rho (R)')
    report = '';
end



function Iq = coreshell(p, q, method)

if method == 1
    if p.sig_a ==0
        a = p.a;
        nr1 = 1;
        dr1 = 1;
    else
        a = linspace(0.01, 1, 15);
        a = a*(p.a + 6*p.a*p.sig_a);
        dr1 = mean(diff(a));
        %nr1 = schultz(p.coreRdivshR, p.sigCore, r1);
        nr1 = schultz(p.a, p.a*p.sig_a, a);
    end
    if p.sig_shell==0
        th = p.shell;
        nr2 = 1;
        dr2 = 1;
    else
        %r2 = linspace(1, p.shellR+3*p.sigshellR, 30);
        th = linspace(0.01, 1, 15);
        th = th*(p.shell + 6*p.shell*p.sig_shell);
        dr2 = mean(diff(th));
        %nr2 = schultz(p.shellR, p.sigshellR, r2);
        nr2 = schultz(p.shell, p.shell*p.sig_shell, th);
    end
    
    yf = zeros(size(q));
    for mm=1:numel(nr2);
        for kk=1:numel(nr1)
            Iq = integ_F_mu(q, a(kk), p.eccentricity, th(mm), ...
                p.core_rho, p.shell_rho, p.solvent_rho);
            yf = yf + Iq*nr2(mm)*nr1(kk);
        end
    end
end

Iq = yf(:)*dr1*dr2;

function F = F_mu(q, b, e, t, rho_c, rho_shell, rho_solv)
% b is the semi-equatorial
% a is semi-principal axis
% b>a : disk-like oblate.
% b<a : needle-like prolate.
u = linspace(0, 1, 30);diffu = diff(u); deltau = mean(diffu);
q = q(:);

a = b*e;
xc = q*b*sqrt(1+u.^2*(e^2-1));
e2 = (a+t)/(b+t);
xt = q*(b+t)*sqrt(1+u.^2*(e2^2-1));
Vc = 4/3*pi*a*b^2;
Vt = 4/3*pi*(a+t)*(b+t)^2;
F = (rho_c - rho_shell)*Vc*3*j1(xc) + (rho_shell - rho_solv)*Vt*3*j1(xt);
F = F*deltau;

function Iq = integ_F_mu(q, b, e, t, rho_c, rho_shell, rho_solv)
A = F_mu(q, b, e, t, rho_c, rho_shell, rho_solv);
Iq = sum(abs(A).^2,2);

function y = j1(x)
y = (sin(x)-x.*cos(x))./x.^3;
