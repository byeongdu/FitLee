function [out, report] = FitLee_schultzsphere_twobackgs(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Schultz polydisperse sphere fit in absolute unit. ' ,...
'I(q) = powI*q.^PorodExp + fn0*(delta_rho0^2*Sq*P(q; r0, sig0) + ',...
'       Nratio1*delta_rho1^2*P1(q; I0_1, s_1, Rg_1, P_1) + ',...
'       Nratio2*delta_rho1^2*P2(q; I0_2, s_2, Rg_2, P_2)) + background',...
'    P(q) = SchultzsphereFun(q, r0, sig0)',...
'    P1(q) = guinierporodmodel(q, I0_1, s_1, Rg_1, P_1)',...
'    P2(q) = guinierporodmodel(q, I0_2, s_2, Rg_2, P_2)',...
'    background = poly1*q.^poly2 + poly3*q + poly4 + SF_userBG1*UsersBackground1 + SF_userBG1*UsersBackground2',...
'Parameters',...
'  fn0 : volume fraction of particle 0 ',...
'  delta_rho0 : electron density difference for particle 0 (A^-3)',...
'  r0 : Schultz size distribution, radius peak (A)',...
'  sig0: Schultz size distribution, FWHM (A)',...
'  Nratio : Number ratio of particle 1 over particle 0',...
'  delta_rho1 : electron density difference for particle 1 (A)',...
'  D : Interparticle distance (A) (hard sphere potential S(q))',...
'  vf : volume fraction of particle 0 (hard sphere potential S(q)) ',...
'  SF_userBG1 : Scale factor for the user input background1',...
'  UsersBackground1 : Two column formatted user"" input background1',...
'  SF_userBG2 : Scale factor for the user input background1',...
'  UsersBackground2 : Two column formatted user"" input background2',...
'  ',...
' P1(q) = guinierporodmodel(q, G, d, Rg, s)',...
'   I0 : Guinier I0',...
'   Rg : Rg of a particle',...
'   P : Porod exponent of the particle (high q)',...
'   s : Dimension of the particle (small q)',... 
'',...
'Note that',...
'I(q)_diff_cross_section = (delta_rho*r_e)^2*f_n*P(q)',...
'         , where P(0) = V_p^2'...
' ',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    1. Kwon et al. Nature Materials, 2015, 14, 215�223. ',...
'    2. Wang et al. J. Phys. Chem. C. 2013. 117(44), 22627. '};

if numel(varargin) > 1
    p = varargin{1};
    q = varargin{2};
    isini = 0;
elseif isscalar(varargin)
    p = varargin{1};
    isini = 1;
    if ishandle(p)
        guiadd_two(p);
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
    %Nf = p;
    bestP = [];
    bestP.delta_rho0 = 1;
    bestP.r0 = 100;
    bestP.sig0 = 10;
    bestP.fn0 = 1;
%    bestP.delta_rho1 = 1;
%    bestP.Nratio1 = 0;
    bestP.I0_1 = 1;
    bestP.s_1 = 0;
    bestP.Rg_1 = 30;
    bestP.P_1 = 1;
%    bestP.Nratio2 = 1;
    % bestP.delta_rho2 = 1;
    bestP.I0_2 = 1;
    bestP.s_2 = 0;
    bestP.Rg_2 = 30;
    bestP.P_2 = 4;
    bestP.D = 100;
    bestP.vf = 0.0;
    bestP.powI = 1E-4;
    bestP.PorodExp = -4;
    % Need 4 parameters for background.
    bestP.poly1 = 0;
    bestP.poly2 = -2;
    bestP.poly3 = 0;
    bestP.poly4 = 0;
    bestP.SF_userBG1 = 1;
    bestP.SF_userBG2 = 1;
    bestP.string = 'FitLee_schultzsphere_twobackgs(figH)';
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
    UBG1 = evalin('base', 'userbackground1');
    UBG1 = interp1(UBG1(:,1), UBG1(:,2), q);
    UBG2 = evalin('base', 'userbackground2');
    UBG2 = interp1(UBG2(:,1), UBG2(:,2), q);
catch
    UBG1 = zeros(size(q));
    UBG2 = zeros(size(q));
end

r_e = 2.818E-5; % Angstrom
Angstrom2Centimeter = 1E-8;
%r_e2 = (r_e*Angstrom2Centimeter)^2;

q = q(:);

[Pq1, V1, V0] = SchultzsphereFun(q, p.r0, p.sig0);
Pq0 = V1*Pq1(:);


Sq = strfactor2(q, p.D, p.vf);
Pq1 = guinierporodmodel(q, 1, p.P_1, p.Rg_1, p.s_1);
Pq2 = guinierporodmodel(q, 1, p.P_2, p.Rg_2, p.s_2);
Pq1 = Pq1*p.I0_1;
Pq2 = Pq2*p.I0_2;

%Sq = p(4)*q(:).^p(5) + strfactor_2Dpara(q, p(6), p(7));
%Iq = p(1)*Imat*nr1/sum(nr1)
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4 + p.SF_userBG1*UBG1 + p.SF_userBG2*UBG2;
pnumberfraction = p.fn0;
f = pnumberfraction*r_e^2;
Sqa = Sq + Pq1;
%Iq0 = f*p.delta_rho0^2*Pq0.*Sq/Angstrom2Centimeter;
Iq0 = f*p.delta_rho0^2*Pq0.*Sqa/Angstrom2Centimeter;
%Iq1 = f*p.Nratio1*p.delta_rho1^2*Pq1/Angstrom2Centimeter;
%Iq2 = f*p.Nratio2*p.delta_rho0^2*Pq2/Angstrom2Centimeter;
Ipw = p.powI./q.^p.PorodExp;
Iq = Iq0;% + Iq1 + Iq2;
out = Ipw + Iq + Pq2 + back; %Iq2;
%out = [out(:), Ipw(:), Iq0(:), Iq1(:), Iq2(:), Ipw(:), back(:)];
out = [out(:), Ipw(:), Iq0(:), Sqa(:), Pq2, back(:)];
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
    fprintf('Volume of particle0 : %0.3e A^3\n', V0);
    pnumberfraction = pnumberfraction/Angstrom2Centimeter^3; % num fraction in 1/cm^3 unit.
    fprintf('Number fraction of particle0 : %0.3e particles/cm^3. \n', pnumberfraction)
    fprintf('Mol concentration : %0.3e M (mole/L). \n', pnumberfraction/6.022E23/1E-3)
    ns = input('Type the  density of a particle in g/mL unit. e.g. water = 1;');
    fprintf('Weight concentration of your particles: %0.3e g/mL. \n', p.fn0*ns)
    
    report = '';
end
function guiadd_two(figH)
    uph = findobj(figH, 'type', 'uipanel');
    pos = get(figH, 'position');
    hFigHeight = pos(end);
    popup = uicontrol(...
          'parent', uph,...
          'horizontalalignment'       , 'Left', ...
          'Style', 'pushbutton',...
          'String', 'Load Background 1',...
          'Position', [50, hFigHeight-235, 150, 50],...
          'Callback', @setshape1);    
    popup = uicontrol(...
          'parent', uph,...
          'horizontalalignment'       , 'Left', ...
          'Style', 'pushbutton',...
          'String', 'Load Background 2',...
          'Position', [50, hFigHeight-295, 150, 50],...
          'Callback', @setshape2);
function setshape1(varargin)
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
    assignin('base', 'userbackground1', data)
function setshape2(varargin)
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
    assignin('base', 'userbackground2', data)