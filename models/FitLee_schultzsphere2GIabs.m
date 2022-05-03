function [out, report] = FitLee_schultzsphere2GIabs(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Schultz polydisperse sphere fit in $d*d\Sigma/d\Omega$. ' ,...
'I(q) = f0*(delta_rho0^2*Sq*P(q; r0, sig0) + Nratio*delta_rho1^2*P(q; r1, sig1)) + background',...
'    Sq = powI*q.^PorodExp + S(q; D, vf)',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'Parameters',...
'  f0 : Number fraction of particle 0 per unit area = N0/S',...
'  delta_rho0 : electron density difference for particle 0 (A^-3)',...
'  r0 : Schultz size distribution, radius peak (A)',...
'  sig0: Schultz size distribution, FWHM (A)',...
'  Nratio : Number ratio of particle 1 over particle 0',...
'  delta_rho1 : electron density difference for particle 1 (A)',...
'  D : Interparticle distance (A) (hard sphere potential S(q))',...
'  vf : volume fraction of particle 0 (hard sphere potential S(q)) ',...
'',...
'Note that',...
'd*I(q)_diff_cross_section = (delta_rho*r_e)^2*N/S*P(q)',...
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
    bestP.f0 = 1;
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

[Pq1, V1, V0] = SchultzsphereFun(q, p.r0, p.sig0);
Pq1 = V1*Pq1(:);
if abs(p.f0) > 0
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
%pnumberfraction = p.f0;
    [~, ~, ~, S] = schultzRg(p.r0, p.sig0);
    S = S/4; % 4 pi R(A) ^2 --> pi*R(nm)^2
f0 = (1-p.f0)*p.f0/S;
out = f0*r_e^2*(p.delta_rho0^2*Pq1.*Sq+ p.Nratio*p.delta_rho1^2*Pq2)+back;
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
    figure;subplot(2,1,1)
    plot(x, nr);xlabel('Radius (A)');ylabel('n(r)')
    subplot(2,1,2)
    Vr = nr(:).*x(:).^3;
    plot(x, Vr);xlabel('Radius (A)');ylabel('V(r)')
    volumedistribution = [x(:), Vr];
    assignin('base', 'volumedistribution', volumedistribution);
    
    fprintf('Volume of particle0 : %0.3e A^3\n', V0);
    %[~, ~, ~, S] = schultzRg(p.r0, p.sig0);
    %S = S/4/100; % 4 pi R(A) ^2 --> pi*R(nm)^2
    pnumberfraction = p.f0 / (S*100);
    fprintf('Number fraction of particle0 : %0.3e/nm^2. \n', pnumberfraction)
    report = '';
end
