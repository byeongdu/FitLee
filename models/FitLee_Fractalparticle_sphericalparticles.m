function [out, report] = FitLee_Fractalparticle_sphericalparticles(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Spherical (Guinier-Porod) primary particles in a fractal cluster.' ,...
'',...
'$I(q) = fn \cdot r_e^2 \cdot S_c(q) \cdot \left( V_{PP}^2 \cdot \delta\rho_{PP}^2 \cdot P_{PP}(q) + N_{ratio,2} \cdot \delta\rho_2^2 \cdot P_2(q) \right) + I_b$',...
'$\qquad    S_c(q) = S_{HS}(q; R_{h,PP}, v_{f,PP}) + P_{cluster}(q)$',...
'$\qquad    P_{cluster}(q) = $ guinierporodmodel$(q, I_{0,c}, P_c, Rg_c, s_c)$',...
'$\qquad    P_{PP}(q) = $ guinierporodmodel$(q, 1, P_{PP}, Rg_{PP}, 0)$ (shape only)',...
'$\qquad    P_2(q) = $ Schultz polydisperse sphere with $V_2^2$ normalization',...
'$\qquad    V_{PP} = (4\pi/3) \cdot (Rg_{PP} \sqrt{5/3})^3$ (sphere-equivalent from $Rg_{PP}$)',...
'$\qquad    I_b = poly1 \cdot q^{poly2} + poly3 \cdot q + poly4 + SF_{userBG} \cdot UBG$',...
'',...
'Note: both populations share the same $S_c(q)$; no separate $S(q)$ for particle 2.',...
'Output is NOT divided by Angstrom2Centimeter (not in true $cm^{-1}$).',...
'',...
'$Parameters$',...
'$\quad  fn$ : overall scale (volume fraction)',...
'$\quad  delta\_rho\_PP$ : electron density contrast for primary particle $(A^{-3})$',...
'$\quad  Rg_{PP}, P_{PP}$ : Guinier-Porod parameters of primary particle shape',...
'$\quad  Rh\_PP, vf\_PP$ : hard-sphere $S(q)$ for PP arrangement',...
'$\quad  I0\_Cluster, s\_c, Rg\_c, P\_c$ : Guinier-Porod parameters of clusters of PP',...
'$\quad  Nratio2$ : number ratio of secondary spherical particles vs PP',...
'$\quad  delta\_rho2$ : electron density contrast for particle 2 $(A^{-3})$',...
'$\quad  r\_sphere, \sigma\_sphere$ : Schultz radius peak (A) and FWHM (A) for particle 2',...
'$\quad  Rh\_sphere, v_{f,sphere}$ : (declared but unused in current code)',...
'$\quad  SF\_userBG$ : scale factor for the user-loaded background',...
'',...
'Note: $r_e$ = 2.818E-5 A is the classical electron radius.',...
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
    bestP.delta_rho_PP = 1;
    bestP.Rg_PP = 300;
    bestP.P_PP = 4;
    bestP.Rh_PP = 100;
    bestP.vf_PP = 0.0;

    % Primary spherical Particles
    bestP.Nratio2 = 0;
    bestP.delta_rho2 = 0.466;
    bestP.r_sphere = 10;
    bestP.sigma_sphere = 1;
    bestP.Rh_sphere = 20;
    bestP.vf_sphere = 0.0;

    % Clusters of the primary particles.
    bestP.I0_cluster = 1;
    bestP.s_cluster = 1;
    bestP.Rg_cluster = 300;
    bestP.P_cluster = 4;
    
    
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
PP = guinierporodmodel(q, 1, p.P_PP, p.Rg_PP, 0);
%PP = guinierporodmodel(q, p.I0_PP, p.P_PP, p.Rg_PP, p.s_PP);

% Spherical Particles
if p.Nratio2 == 0 || p.r_sphere == 0
    Pq2 = 1;
    Sq2 = 1;
    V0 = 1;
else
    [Pq2, V2, V0] = SchultzsphereFun(q, p.r_sphere, p.sigma_sphere);
    Pq2 = V2*Pq2(:);
    Sq2 = strfactor2(q, p.Rh_sphere, p.vf_sphere);
end

% Fractal Primary Particles
volfraction = p.fn;
f = volfraction*r_e^2;
% Rg = sqrt(3/5)*R. Therefore, R = Rg*sqrt(5/3) and volume = 4*pi/3*R^3
R = p.Rg_PP*sqrt(5/3);
PPvol = 4*pi/3*R^3;
%Iq1 = f*PPvol^2*PP.*Sq;
Iq1 = PPvol^2*(p.delta_rho_PP^2*PP);
%Iq1 = f*PPvol^2*(p.delta_rho_PP^2*PP + Pcluster);
% Spherical Particles
Iq2 = p.Nratio2*p.delta_rho2^2*Pq2.*Sq2;

% background
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4 + p.SF_userBG*UBG;


% Intensity
out = f*Sq.*(Iq1 + Iq2) + back;
if numel(Iq2) == 1
    Iq2 = Iq2*ones(size(Iq1));
end
out = [out(:), Iq1(:), Iq2(:), back(:)];
%out = out*Angstrom2Centimeter; % data is in cm^-1
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