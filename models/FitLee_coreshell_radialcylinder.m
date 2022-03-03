function [out, report] = FitLee_coreshell_radialcylinder(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Schultz polydisperse coreshell cylinder radial direction fit. ' ,...
    'Data needs to be plotted as q vs Iq*q',...
'I(q) = I0*P(q; r0, sig0)*S(q) + background',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
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
    bestP.core_rho = 0;
    bestP.shell_rho = 0.39;
    bestP.solvent_rho = 0.3344;
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
        edensity = [p.core_rho, p.shell_rho, p.solvent_rho];
        Iq0 = saxscylinder_CS(q, radius, edensity, 0);
        S = S + (pi*R)^2*nr(i)*nr2(j);
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

pnumberfraction = p.volf/(Vcyl*Angstrom2Centimeter^3);
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;

out = pnumberfraction*r_e^2*Angstrom2Centimeter^2*Iq.*Sq + back;
if isnan(out)
    out = ones(size(out))*1E20;
end

if nargout == 2
    maxR = p.coreR+p.shellthick;
    x = linspace(maxR/10, maxR + maxR*sigma*9, 100);
    nr = schultzdist(x, maxR, maxR*p.sigma);

    
    figure;
    plot(x, nr);xlabel('Radius (A)');ylabel('n(R)')
    report = '';
end


