function [out, report] = FitLee_helicalcylinder(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Schultz polydisperse helical cylinder (spring) radial direction fit. ' ,...
    'Data needs to be plotted as q_r vs I(q)',...
'I(q) = \phi(\Delta{\rho})^2P(q_r; R, R_h, p_h, Npitch or H) + I_b',...
'    \phi or volf : volume fraction of the spring',...
'    \Delta{\rho} : Contrast or electron density difference (#e/A^3)',...
'    P(q_r) : form factor and P(0) = Vp^2',...
'    R : Radius of cylinder',...
'    R_h : Radius of helix (0<=R_h<=\surd(2\pi)p_h)',...
'    p_h : Pitch length',...
'    Npitch : Number of the pitch',...
'    H = Npitch \cdot p_h is the spring length',...
'    Use Npitch or H. Then set the other 0',...
'    Npitch \cdot sqrt((2\pi R_h)^2+p_h^2) is the contour length of the spring',...
'    I_b = poly1*q.^{poly2} + poly3*q + poly4',...
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
    bestP.delta_rho = 0;
    bestP.R = 25;
    bestP.R_frac = 0.1;
    bestP.Rh = 20;
    bestP.Rh_frac = 0;
    bestP.ph = 100;
    bestP.Npitch = 1;
    bestP.H = 0;

    bestP.poly1 = 0;
    bestP.poly2 = 0;
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

% radial direction
q = q(:);
numberofPoint = 10;
if p.R_frac > 0
    [nr, x] = schultzdist99(p.R, p.R_frac*p.R, numberofPoint);
else
    nr = 1;
    x = p.R;
end

if p.Rh_frac > 0
    [nrRh, xRh] = schultzdist99(p.Rh, p.Rh_frac*p.Rh, 10);
else
    nrRh = 1;
    xRh = p.Rh;
end

V = 0;

if p.H==0
    N = p.Npitch;
else
    N = p.H/p.ph;
end

v = p.ph/(2*pi*p.Rh);
vnz = (1+v^2)^(-1/2);
Rh_vnz = p.Rh/vnz;
helixlength = 2*pi*Rh_vnz;
%vol = helixlength*pi*p.R^2;

Iq = zeros(size(q));
Iqa = Iq;
for k=1:numel(xRh)
    R_h = xRh(k);
    vnz = 2*pi*R_h/helixlength;
    v = sqrt(vnz^-2-1);
    ph = 2*pi*R_h*v;

    for i=1:numel(x)
        R = x(i);
        [FF, Vp] = saxs_spiralcylinder(q, zeros(size(q)), zeros(size(q)), [R_h, R, ph, N]);
        Iq = Iq + nr(i)*abs(FF).^2;
        V = V + Vp*nr(i);
    end
    Iqa = Iqa + Iq*nrRh(k);
end


back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;

pnumberfraction = p.volf/(V*Angstrom2Centimeter^3);
out = [pnumberfraction*r_e^2*Angstrom2Centimeter^2*p.delta_rho^2*Iqa+back, back];
if isnan(out)
    out = ones(size(out));
end



if nargout == 2
    maxR = p.R;
    x = linspace(maxR/10, maxR + maxR*p.R_frac*9, 100);
    nr = schultzdist(x, maxR, maxR*p.R_frac);

    
    figure;
    plot(x, nr);xlabel('Radius (A)');ylabel('n(R)')
    report = '';
end


