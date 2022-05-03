function [out, report] = FitLee_3level(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
global FitLee_helpstr
FitLee_helpstr = {'Schultz polydisperse sphere fit in absolute unit. ' ,...
'I(q) = powI*q.^PorodExp + fn0*delta_rho0^2*P0(q; r0, sig0)*Sq0 + background',...
'    Sq0 = Sq1*P2(q) + S_inf(q; D, vf)',... % structure factor of primary particle
'    P1(q) = guinierporodmodel(q, I0_1, s_1, Rg_1, P_1)',... % form factor of level 2 structure
'    S2q = (I3(q) + S2_inf(q))',... % structure factor of level 2 structure.
'    S2_inf(q) = 1',... % Ordering of level 2 structure.
'    I3(q) = Nratio*delta_rho1^2*P(q; r1, sig1) + UnifiedSAXS(q; G, B, Rg, P)',...
'    background = poly1*q.^poly2 + poly3*q + poly4 + SF_userBG*UsersBackground',...
'Parameters',...
'  fn0 : volume fraction of particle 0 ',...
'  delta_rho0 : electron density difference for particle 0 (A^-3)',...
'  r0 : Schultz size distribution, radius peak (A)',...
'  sig0: Schultz size distribution, FWHM (A)',...
'  Nratio : Number ratio of particle 1 over particle 0',...
'  delta_rho1 : electron density difference for particle 1 (A)',...
'  D : Interparticle distance (A) (hard sphere potential S(q))',...
'  vf : volume fraction of particle 0 (hard sphere potential S(q)) ',...
'  SF_userBG : Scale factor for the user input background',...
'  UsersBackground : Two column formatted user"" input background',...
'  ',...
' P1(q) = guinierporodmodel(q, G, d, Rg, s)',...
'   I0 : Guinier I0',...
'   Rg : Rg of a particle',...
'   P : Porod exponent',...
'   s : dimension of the particle',... 
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
    if ishandle(p)
        guiadd(p);
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
    bestP.G1 = 10000;
    bestP.B1 = 1E-6;
    bestP.Rg1 = 300;
    bestP.P1 = 4;
    bestP.I0_2 = 1;
    bestP.s_2 = 2;
    bestP.Rg_2 = 300;
    bestP.P_2 = 4;
    bestP.D = 100;
    bestP.vf = 0.0;
    bestP.powI = 0;
    bestP.PorodExp = -4;
    % Need 4 parameters for background.
    bestP.poly1 = 0;
    bestP.poly2 = -2;
    bestP.poly3 = 0;
    bestP.poly4 = 0;
    bestP.SF_userBG = 1;
    bestP.string = 'FitLee_3level(figH)';
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


try
    UBG = evalin('base', 'userbackground');
    UBG = interp1(UBG(:,1), UBG(:,2), q);
catch
    UBG = zeros(size(q));
end

r_e = 2.818E-5; % Angstrom
Angstrom2Centimeter = 1E-8;
%r_e2 = (r_e*Angstrom2Centimeter)^2;

q = q(:);

[Pq1, V1, V0] = SchultzsphereFun(q, p.r0, p.sig0);
Pq1 = V1*Pq1(:);
if abs(p.fn0) > 0
    [Pq2, V2] = SchultzsphereFun(q, p.r1, p.sig1);
    Pq2 = V2*Pq2(:);
else
    Pq2 = 0;
end

Sq1 = strfactor2(q, p.D, p.vf);
Pq_1 = guinierporodmodel(q, p.I0_2, p.P_2, p.Rg_2, p.s_2);
Pq_unfied = p.G1*(exp(-q.^2*p.Rg1^2/3) + p.B1*(erf(q*p.Rg1/sqrt(6)).^3./q).^p.P1);

Poq = p.powI*q.^p.PorodExp;

Sq = (Pq_unfied + Poq + 1).*Pq_1 + Sq1;
%Sq = p(4)*q(:).^p(5) + strfactor_2Dpara(q, p(6), p(7));
%Iq = p(1)*Imat*nr1/sum(nr1)
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4 + p.SF_userBG*UBG;
pnumberfraction = p.fn0/V0;
%out = Poq + pnumberfraction*r_e^2*Angstrom2Centimeter^2*(p.delta_rho0^2*Pq1.*Sq+ p.Nratio*p.delta_rho1^2*Pq2)+back;
out = pnumberfraction*r_e^2*(p.delta_rho0^2*Pq1.*Sq + p.Nratio*p.delta_rho1^2*Pq2)/Angstrom2Centimeter+back;
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
    pnumberfraction = pnumberfraction/Angstrom2Centimeter^3; % num fraction in 1/cm^3 unit.
    fprintf('Volume of particle0 : %0.3e A^3\n', V0);
    fprintf('Number fraction of particle0 : %0.3e particles/cm^3. \n', pnumberfraction)
    fprintf('Mol concentration : %0.3e mol/mL. \n', pnumberfraction/6.022E23)
    ns = input('Type the  density of a particle in g/mL unit. e.g. water = 1;');
    fprintf('Weight concentration of your particles: %0.3e g/mL. \n', p.fn0*ns)
    
    report = '';
end
function guiadd(figH)
    uph = findobj(figH, 'type', 'uipanel');
    pos = get(figH, 'position');
    hFigHeight = pos(end);
        uicontrol(...
          'style'                     , 'text', ...
          'backgroundcolor'           , [0,0.5,0.5], ...
          'foregroundcolor'           , [1,1,1], ...
          'horizontalalignment'       , 'Left', ...
          'parent'                    , uph, ...
          'string'                  ,   'BackG',...
          'position'                  ,[1, hFigHeight-210,100,20], ...
          'tag'                       , 'text_Pshape');
    popup = uicontrol(...
          'parent', uph,...
          'horizontalalignment'       , 'Left', ...
          'Style', 'pushbutton',...
          'String', 'Load Background',...
          'Position', [50, hFigHeight-235, 150, 50],...
          'Callback', @setshape);
function setshape(varargin)
    obj = varargin{1};
    [filename, pathname] = uigetfile('*.txt; *.dat; *.*', 'Pick your file');
    if isequal(filename,0) || isequal(pathname,0)
        disp('User pressed cancel')
        return
    else
        disp(['User selected ', fullfile(pathname, filename)])
    end
    filename = fullfile(pathname, filesep, filename);
    [~, data] = hdrload(filename);
    assignin('base', 'userbackground', data)