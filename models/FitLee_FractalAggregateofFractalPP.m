function [out, report] = FitLee_FractalAggregateofFractalPP(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Fractal aggregates of fractal primary particles, $(cm^{-1})$.' ,...
'',...
'$I(q) = (fn \cdot r_e^2 / 10^{-8}) \cdot \left( V_{PP}^2 \cdot P_{PP}(q) \cdot S_c(q) + N_{ratio,2} \cdot \delta\rho_2^2 \cdot P_2(q) \cdot S_2(q) \right) + I_b$',...
'$\qquad    S_c(q) = S_{HS}(q; R_{h,PP}, v_{f,PP}) + P_{cluster}(q)$',...
'$\qquad    P_{cluster}(q) = $ guinierporodmodel$(q, I_{0,c}, P_c, Rg_c, s_c)$',...
'$\qquad    P_{PP}(q) = $ guinierporodmodel$(q, I_{0,PP}, P_{PP}, Rg_{PP}, 0)$',...
'$\qquad    P_2(q) = $ Schultz polydisperse sphere with $V_2^2$ normalization',...
'$\qquad    S_2(q) = S_{HS}(q; R_{h,2}, v_{f,2})$',...
'$\qquad    V_{PP} = (4\pi/3) \cdot (Rg_{PP} \sqrt{5/3})^3$ (sphere-equivalent from $Rg_{PP}$)',...
'$\qquad    I_b = poly1 \cdot q^{poly2} + poly3 \cdot q + poly4 + SF_{userBG} \cdot UBG$',...
'',...
'Note: $I_{0,PP}$ (inside $P_{PP}$) absorbs the role of $\delta\rho_{PP}^2$ for the primary particle.',...
'$P_2$ uses its own structure factor $S_2$; both populations share the same $fn$ scale.',...
'',...
'$Parameters$',...
'$\quad  fn$ : overall scale (volume fraction)',...
'$\quad  I0\_PP, P_{PP}, Rg_{PP}$ : Guinier-Porod parameters of primary particle (PP)',...
'$\quad  Rh\_PP, vf\_PP$ : hard-sphere $S(q)$ for PP arrangement',...
'$\quad  I0\_Cluster, s\_c, Rg\_c, P\_c$ : Guinier-Porod parameters of clusters of PP',...
'$\quad  Nratio2$ : number ratio of secondary spherical particles vs PP',...
'$\quad  \delta\rho\_2$ : electron density contrast for particle 2 $(A^{-3})$',...
'$\quad  r2, \sigma\_2$ : Schultz radius peak (A) and FWHM (A) for particle 2',...
'$\quad  Rh2, vf2$ : hard-sphere $S(q)$ for particle 2 arrangement',...
'$\quad  SF\_userBG$ : scale factor for the user-loaded background',...
'',...
'Guinier-Porod sub-parameter meanings:',...
'$\quad  I\_0$ : Guinier prefactor',...
'$\quad  Rg$ : radius of gyration (A)',...
'$\quad  P$ : Porod exponent at high q',...
'$\quad  s$ : dimensionality at low q (0 sphere, 1 rod, 2 plate)',...
'',...
'Note: $r_e$ = 2.818E-5 A is the classical electron radius.',...
'$  d\Sigma/d\Omega = (\delta\rho \cdot r_e)^2 \cdot f_n \cdot P(q),$ where $P(0) = V_p^2$',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
};

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
    bestP.fn = 1;
    % Primary Fractal Particles
    bestP.I0_PP = 1;
    bestP.Rg_PP = 300;
    bestP.P_PP = 4;
    % Clusters of the primary particles.
    bestP.Rh_PP = 100;
    bestP.vf_PP = 0.0;
    bestP.I0_cluster = 1;
    bestP.s_cluster = 1;
    bestP.Rg_cluster = 300;
    bestP.P_cluster = 4;
    
    % Secondary Particles
    bestP.Nratio2 = 0;
    bestP.delta_rho2 = 0.466;
    bestP.r2 = 10;
    bestP.sigma2 = 1;
    bestP.Rh2 = 20;
    bestP.vf2 = 0.0;
    
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

% Structure factor for Primary Particles
Sq0 = strfactor2(q, p.Rh_PP, p.vf_PP);
Pcluster = guinierporodmodel(q, p.I0_cluster, p.P_cluster, p.Rg_cluster, p.s_cluster);
Sq = Sq0 + Pcluster;
% Primary Particles
PP = guinierporodmodel(q, p.I0_PP, p.P_PP, p.Rg_PP, 0);
%PP = guinierporodmodel(q, p.I0_PP, p.P_PP, p.Rg_PP, p.s_PP);


% Spherical Particles
if p.Nratio2 == 0 || p.r2 == 0
    Pq2 = 1;
    Sq2 = 1;
    V0 = 1;
else
    [Pq2, V2, V0] = SchultzsphereFun(q, p.r2, p.sigma2);
    Pq2 = V2*Pq2(:);
    Sq2 = strfactor2(q, p.Rh2, p.vf2);
end

% Fractal Primary Particles
pnumberfraction = p.fn;
f = pnumberfraction*r_e^2/Angstrom2Centimeter;
% Rg = sqrt(3/5)*R. Therefore, R = Rg*sqrt(5/3) and volume = 4*pi/3*R^3
R = p.Rg_PP*sqrt(5/3);
PPvol = 4*pi/3*R^3;
Iq1 = f*PPvol^2*PP.*Sq;
%Iq1 = f*PPvol^2*(PP + Pcluster);
% Spherical Particles
Iq2 = f*p.Nratio2*p.delta_rho2^2*Pq2.*Sq2;

% background
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4 + p.SF_userBG*UBG;


% Intensity
out = Iq1 + Iq2 + back;
if numel(Iq2) == 1
    Iq2 = Iq2*ones(size(Iq1));
end
out = [out(:), Iq1(:), Iq2(:), back(:)];
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