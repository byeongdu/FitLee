function [out, report] = FitLee_polyhedra(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Schultz polydisperse polyhedra fit. ' ,...
'I(q) = I0*(Sq*P(q; edgelength, sig0) + I1*P(q; r1, sig1)) + background',...
'    Sq = powI*q.^PorodExp + S(q; D, vf)',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
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
    if ishandle(p)
        guiadd(p);
        return
    end
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
    bestP.fn0 = 1;
    bestP.delta_rho0 = 4.66;
    bestP.edgelength = 100;
    bestP.sig0 = 10;
    bestP.Nratio = 0;
    bestP.delta_rho1 = 4.66;
    bestP.r1 = 0;
    bestP.sig1 = 0;
    bestP.D = 100;
    bestP.vf = 0.0;
    bestP.powI = 0;
    bestP.PorodExp = -4;
    % Need 4 parameters for background.
    bestP.poly1 = 0;
    bestP.poly2 = 0;
    bestP.poly3 = 0;
    bestP.poly4 = 0;
    bestP.string = 'FitLee_polyhedra(figH)';
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

% special for cube particle .==========
%dttemp = evalin('base', 'cube120');
try
    fs = evalin('base', 'formfactortoscale');
    dttemp = fs.data;
    dtsize = fs.size;
catch
    sprintf('The format of formfactortoscale is:')
    sprintf('           formfactortoscale.data')
    sprintf('           formfactortoscale.size')
    error('formfactortoscale is not defined')
end

if iscell(fs.size)
    r2 = fs.size{2};
    dtsize = fs.size{1};
    nr2 = schultzdist(r2, p.r1, p.sig1); % distribution of the other dimension.
    dt = sum(dttemp(:,2:end).*repmat(nr2(:)', numel(dttemp(:,1)), 1)*(r2(2)-r2(1)), 2);
    dtq = dttemp(:,1);
    dttemp = [dtq, dt];
end


Pq1 = formfactorscale(dttemp(:,1), dttemp(:,2), dtsize, q, p.edgelength, p.sig0, 1);
Pq1 = Pq1(:);
V0 = sqrt(dttemp(1, 2));

if abs(p.Nratio) > 0
    [Pq2, V2] = SchultzsphereFun(q, p.r1, p.sig1);
    Pq2 = V2*Pq2(:);
else
    Pq2 = 0;
end

pnumberfraction = p.fn0/(V0);
f = pnumberfraction/Angstrom2Centimeter;

Sq1 = strfactor2(q, p.D, p.vf);
Sq = p.powI*q.^p.PorodExp + Sq1;
%Sq = p(4)*q(:).^p(5) + strfactor_2Dpara(q, p(6), p(7));
%Iq = p(1)*Imat*nr1/sum(nr1)
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
out = f*r_e^2*(p.delta_rho0^2*Pq1.*Sq+ p.Nratio*p.delta_rho1^2*Pq2)+back;
%out = pnumberfraction*r_e^2*Angstrom2Centimeter^2*(p.delta_rho0^2*Pq1.*Sq+ p.Nratio*p.delta_rho1^2*Pq2)+back;

if isnan(out)
    out = ones(size(out))*1E20;
end

if nargout == 2
    sig = max(p.sig0, p.sig1)*10;
    r0 = max(p.edgelength, p.r1);
    if sig < r0*2
        sig = r0;
    end
    x = 0:1:(r0+sig);
    nr0 = schultz(p.edgelength, p.sig0, x);
    nr = nr0;
    if p.Nratio > 0
        nr1 = schultz(p.r1, p.sig1, x);
        nr = nr+p.Nratio*nr1;
    end
    figure;
    plot(x, nr);xlabel('Edgelength (A)');ylabel('n(r)')
    fprintf('Volume of particle0 : %0.3e A^3\n', V0);
    fprintf('Number fraction of particle0 : %0.3e particles/cm^3. \n', pnumberfraction)
    report = '';
end

function guiadd(figH)
    uph = findobj(figH, 'type', 'uipanel');
    pos = get(figH, 'position');
    qmax_position = get(findobj(figH, 'string', 'q max'), 'position');
    qmax_position = qmax_position(2);
    hFigHeight = pos(end);
        uicontrol(...
          'style'                     , 'text', ...
          'backgroundcolor'           , [0,0.5,0.5], ...
          'foregroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'Left', ...
          'parent'                    , uph, ...
          'string'                  ,   'Particle Shape :',...
          'position'                  ,[1, qmax_position-30,150,20], ...
          'tag'                       , 'text_Pshape');
    popup = uicontrol(...
          'parent', uph,...
          'horizontalalignment'       , 'Left', ...
          'Style', 'popup',...
          'String', {'cube','octahedron','rhombicdodecahedron','cubooctahedron', 'truncatedcube', 'truncatedoctahedron', 'icosahedron'},...
          'Position', [100, qmax_position-30, 100, 20],...
          'Callback', @setshape);
function setshape(varargin)
    obj = varargin{1};
    list = get(obj, 'string');
    switch list{get(obj, 'value')}
        case 'cube'
            load cube.mat;
        case 'octahedron'
            load octahedron.mat;
        case 'rhombicdodecahedron'
            load rhombicdodecahedron.mat;
        case 'cubooctahedron'
            load cubooctahedron.mat;
        case 'truncatedcube'
            load truncatedcube.mat;
        case 'truncatedoctahedron'
            load truncatedoctahedron.mat;
        case 'icosahedron'
            load icosahedron.mat;
        case 'longRD'
            load longRD.mat;
        case 'convexcube_0.2154'
            load convexcube.mat;
        case 'concavecube_0.1538'
            load concavecube.mat;
        case 'concavecube_0.0769'
            load concavecube2.mat;
        case 'concavecube_0.1'
            load concavecube3.mat;
        case 'concavecube_0.05'
            load concavecube4.mat;
    end
    assignin('base', 'formfactortoscale', formfactortoscale) 

    
