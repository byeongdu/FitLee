function [out, report] = FitLee_truncatedcube(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Schultz polydisperse polyhedra fit. ' ,...
'I(q) = fn0*delta_rho^2*Sq*\int_0^inf n(x; edgelength, sig)*P(q; x, gamma)dx + background',...
'    Sq = powI*q.^PorodExp + S(q; D, vf)',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'Parameters -',...
' edgelength : edgelength of a mother cube.',...
' sig : standard deviation of the cube edgelength.',...
' gamma : fraction of the edgelength to be cut out.'...
'         gamma*edgelength*sqrt(2) is the base length of the truncated tetrahedron.',...
' n(x; edgelength, sig) : Schultz distribution function.',...
'    its peak at edgelength.',...
'    its breadth is sig.',...
' S(q) is the hard sphere structure factor.',...
'     D : hard sphere radius.',...
'     vf : hard sphere volume fraction.',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    1. A Senesi and B Lee, J. Appl. Cryst. 2015. '};
    
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
    bestP.fn0 = 1E-9;
    bestP.delta_rho = 4.66;
    bestP.edgelength = 300;
    bestP.gamma = 0.25;
    bestP.sig = 10;
    bestP.D = 300;
    bestP.vf = 0.0;
    bestP.powI = 0;
    bestP.PorodExp = -4;
    % Need 4 parameters for background.
    bestP.poly1 = 0;
    bestP.poly2 = -2;
    bestP.poly3 = 0;
    bestP.poly4 = 0;
    bestP.string = 'FitLee_truncatedcube(figH)';
    assignin('base', 'bestP', bestP);
    out = bestP;
    return
end

%% fitting codes .......................
if iscell(q)
    if numel(q) > 1
        error('FitLee_polyhedra.m is for fitting a set of data, for now')
    end
    q = q{1};
end
q = q(:);

r_e = 2.818E-5; % Angstrom
Angstrom2Centimeter = 1E-8;

if (p.gamma > 0.5)
    error('gamma should not be more than 0.5.')
end

if (p.gamma < 0)
    error('gamma should not be positive.')
end

% special for cube particle .==========
%dttemp = evalin('base', 'cube120');

if p.sig>0
    numberofPoint = 11;
    [nr, x] = schultzdist99(p.edgelength, p.sig, numberofPoint);
    x = x/p.edgelength;
    Pq = saxs_average(q, 'saxstruncatedcube', x'*[p.edgelength, p.gamma*p.edgelength*sqrt(2)], nr(:), 1);
else
    Pq = saxs_average(q, 'saxstruncatedcube', [p.edgelength, p.gamma*p.edgelength*sqrt(2)]);
end


Sq1 = strfactor2(q, p.D, p.vf);
Sq = p.powI*q.^p.PorodExp + Sq1;

back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
f = p.fn0*r_e^2/Angstrom2Centimeter;
out = f*p.delta_rho^2*Pq.*Sq+back;
%out = pnumberfraction*r_e^2*Angstrom2Centimeter^2*(p.delta_rho0^2*Pq1.*Sq+ p.Nratio*p.delta_rho1^2*Pq2)+back;

if isnan(out)
    out = ones(size(out))*1E20;
end

if nargout == 2
    x = 0:1:(p.edgelength+p.sig)*10;
    nr0 = schultz(p.edgelength, p.sig, x);
    nr = nr0;
    figure;
    plot(x, nr);xlabel('Edgelength (A)');ylabel('n(r)')
    fprintf('Volume of the TC : %0.3e A^3\n', V0);
    fprintf('Number fraction of TC : %0.3e particles/cm^3. \n', pnumberfraction)
    report = '';
end
