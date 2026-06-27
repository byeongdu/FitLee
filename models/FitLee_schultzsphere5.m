function [out, report] = FitLee_schultzsphere5(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Schultz sphere + two Guinier-Porod components, absolute unit $(cm^{-1}). $' ,...
'',...
'$I(q) = I_{pw} + (fn0/\left<V_0\right>) \cdot r_e^2 \cdot \left( \delta\rho_0^2 \cdot P_0(q) \cdot S(q) + N_{ratio,1} \cdot \delta\rho_1^2 \cdot P_1(q) + N_{ratio,2} \cdot \delta\rho_2^2 \cdot P_2(q) \right) + I_b$',...
'$\qquad    I_{pw} = powI \cdot q^{PorodExp}$',...
'$\qquad    P_0(q) = $ SchultzsphereFun$(q, r_0, sig_0)$',...
'$\qquad    P_1(q) = $ guinierporodmodel$(q, I_{0,1}, P_1, Rg_1, s_1)$',...
'$\qquad    P_2(q) = $ guinierporodmodel$(q, I_{0,2}, P_2, Rg_2, s_2)$',...
'$\qquad    S(q) = S_{HS}(q; D, v_f)$',...
'$\qquad    I_b = poly1 \cdot q^{poly2} + poly3 \cdot q + poly4 + SF_{userBG} \cdot UBG$',...
'',...
'Note: $S(q)$ multiplies only $P_0$; the two Guinier-Porod populations are added bare.',...
'All three populations share the same $\left<V_0\right>$ denominator from particle 0.',...
'',...
'$Parameters$',...
'$\quad  fn0$ : volume fraction of particle 0 (dimensionless)',...
'$\quad  delta\_rho0, delta\_rho1, delta\_rho2$ : electron density contrasts $(A^{-3})$',...
'$\quad  r0, sig0$ : Schultz radius peak (A) and FWHM (A) for particle 0',...
'$\quad  Nratio1, Nratio2$ : number ratios for populations 1 and 2 over particle 0',...
'$\quad  I0\_1, P\_1, Rg\_1, s\_1$ : Guinier-Porod parameters for population 1',...
'$\quad  I0\_2, P\_2, Rg\_2, s\_2$ : Guinier-Porod parameters for population 2',...
'$\quad  D$ : interparticle distance (A) for hard-sphere $S_{HS}(q)$',...
'$\quad  vf$ : volume fraction in hard-sphere $S_{HS}(q)$',...
'$\quad  powI, PorodExp$ : standalone low-q power-law term (outside the bracket)',...
'$\quad  SF\_userBG$ : scale factor for the user-loaded background',...
'',...
'Guinier-Porod sub-parameter meanings:',...
'$\quad  I\_0$ : Guinier $I\_0$',...
'$\quad  Rg$ : radius of gyration (A)',...
'$\quad  P$ : Porod exponent at high q',...
'$\quad  s$ : dimensionality at low q (0 sphere, 1 rod, 2 plate)',...
'',...
'Note: $r_e$ = 2.818E-5 A is the classical electron radius.',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    1. Kwon et al. Nature Materials, 2015, 14, 215-223. ',...
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
    bestP.delta_rho1 = 1;
    bestP.Nratio1 = 0;
    bestP.I0_1 = 1;
    bestP.s_1 = 1;
    bestP.Rg_1 = 300;
    bestP.P_1 = 4;
    bestP.Nratio2 = 0;
    bestP.delta_rho2 = 1;
    bestP.I0_2 = 1;
    bestP.s_2 = 1;
    bestP.Rg_2 = 300;
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
    bestP.SF_userBG = 1;
    bestP.string = 'FitLee_schultzsphere5(figH)';
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
Pq0 = V1*Pq1(:);

Sq = strfactor2(q, p.D, p.vf);
Pq1 = guinierporodmodel(q, p.I0_1, p.P_1, p.Rg_1, p.s_1);
Pq2 = guinierporodmodel(q, p.I0_2, p.P_2, p.Rg_2, p.s_2);

%Sq = p(4)*q(:).^p(5) + strfactor_2Dpara(q, p(6), p(7));
%Iq = p(1)*Imat*nr1/sum(nr1)
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4 + p.SF_userBG*UBG;
pnumberfraction = p.fn0/(V0);
f = pnumberfraction*r_e^2;
Iq0 = f*p.delta_rho0^2*Pq0.*Sq/Angstrom2Centimeter;
Iq1 = f*p.Nratio1*p.delta_rho1^2*Pq1/Angstrom2Centimeter;
Iq2 = f*p.Nratio2*p.delta_rho2^2*Pq2/Angstrom2Centimeter;
Ipw = p.powI*q.^p.PorodExp;
Iq = Iq0 + Iq1 + Iq2;
out = Ipw + Iq + back;
out = [out(:), Ipw(:), Iq0(:), Iq1(:), Iq2(:), Ipw(:), back(:)];
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