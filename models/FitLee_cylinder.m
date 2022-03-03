function [out, report] = FitLee_cylinder(varargin)
    % schultz polydisperse radial direction cylinder with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Schultz polydisperse cylinder radial direction fit. ' ,...
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
    bestP.volfrac = 1;
    bestP.delta_rho = 0.322;
    bestP.R = 25;
    bestP.R_frac = 0.1;
    bestP.H = 100;
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

if p.R_frac > 0
    numpnt = 41;
    [nr, x] = schultzdist99(p.R, p.R*p.R_frac, numpnt);
else
    numpnt = 1;
    nr = 1;
    x = p.R;
end
Iq = zeros(size(q));
for i=1:numpnt
    [F, Vcyl] = saxscylinder3(q, x(i));
    Iq = Iq + nr(i)*abs(F).^2*p.H^2;
end
Vcyl = Vcyl * p.H;

r_e = 2.818E-5; % Angstrom
Angstrom2Centimeter = 1E-8;

Sq = strfactor2(q, p.D, p.vf);

pnumberfraction = p.volfrac/(Vcyl*Angstrom2Centimeter^3);
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;

out = pnumberfraction*r_e^2*Angstrom2Centimeter^2*p.delta_rho^2*Iq.*Sq + back;


if isnan(out)
    out = ones(size(out))*1E20;
end

if nargout == 2
    maxR = p.R;
    x = linspace(maxR/10, maxR + maxR*sigma*9, 100);
    nr = schultzdist(x, maxR, maxR*p.sigma);

    
    figure;
    plot(x, nr);xlabel('Radius (A)');ylabel('n(R)')
    report = '';
end


