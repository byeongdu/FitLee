function [out, report] = FitLee_core_fuzzycorona_sphere(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Core-corona micelle with fuzzy corona, absolute unit $(cm^{-1}). $' ,...
'',...
'$I(q) = (fn0/\left<V_m\right>) \cdot r_e^2 \cdot \left( P_{cs}(q) \cdot S(q; D, v_f) + P_c(q) \right) + I_b$',...
'$\qquad    P_{cs}(q) = \left| (\rho_{sh} - \rho_{solv}) V_o P_s(q R_o) e^{-(q \sigma)^2/2} - (\rho_{sh} - \rho_{core}) V_c P_s(q R_o v) \right|^2$',...
'$\qquad    P_c(q) = N_c \cdot (\delta\rho_{chain} \cdot V_{ch})^2 \cdot P_{chain}(q R_g)$',...
'$\qquad    v = $ coreRdivshR = $R_c / R_o$, $\sigma = $ sig\_fuzzy (Gaussian gradient),',...
'$\qquad    P_s(\cdot)$ is the sphere shape factor, $P_{chain}$ is the Debye chain function',...
'$\qquad    I_b = poly1 \cdot q^{poly2} + poly3 \cdot q + poly4 + SF_{userBG} \cdot UBG$',...
'',...
'Note: $S(q; D, v_f)$ multiplies the core-shell term only; the chain term is added bare.',...
'Both terms share the same $\left<V_m\right>$ denominator (from the coreshell helper).',...
'',...
'$Parameters$',...
'$\quad  fn0$ : volume fraction of particles (dimensionless)',...
'$\quad  core\_rho, sh\_rho, solvent\_rho$ : electron densities $(A^{-3})$',...
'$\quad  delta\_rho\_chain$ : excess electron density of chain $(A^{-3})$',...
'$\quad  shellR, sigshellR$ : outer radius $R\_o$ (A) and Schultz relative width',...
'$\quad  coreRdivshR, sigCore$ : core/shell ratio $v$ and its relative width',...
'$\quad  sig\_fuzzy$ : Gaussian density-gradient width $\sigma$ (A) at shell/solvent interface',...
'$\quad  Nc$ : number of chains per particle',...
'$\quad  Rg$ : chain radius of gyration (A)',...
'$\quad  D, vf$ : hard-sphere $S(q)$ parameters',...
'$\quad  SF\_userBG$ : scale factor for the user-loaded background',...
'',...
'Note: $r_e$ = 2.818E-5 A is the classical electron radius.',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    Pedersen and Gerstenberg, J. Chem. Phys., http://aip.scitation.org/doi/pdf/10.1063/1.1493771',...
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
    bestP.fn0 = 1;
    bestP.core_rho = 0.3;
    bestP.coreRdivshR = 0.5;
    bestP.sigCore = 0.1;
    
    bestP.sh_rho = 0.4;
    bestP.shellR = 50;
    bestP.sigshellR = 0.1;
    bestP.sig_fuzzy = 10;
    
    bestP.solvent_rho = 0.3344;
    
    bestP.Nc = 50;
    bestP.delta_rho_chain = 0.1;
    bestP.Rg = 20;
    
    bestP.D = 100;
    bestP.vf = 0.1;

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
[P, Pc, Vm] = coreshell(p, q);

Sq = strfactor2(q, p.D, p.vf);

back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4 + p.SF_userBG*UBG;
pnumberfraction = p.fn0/(Vm);
f = pnumberfraction*r_e^2/Angstrom2Centimeter;
Iq0 = f*P.*Sq;
Iq1 = f*Pc;
out = Iq0 + Iq1 + back;
out = [out(:), Iq0(:), Iq1(:), back(:)];
if isnan(out)
    out = ones(size(out));
end

function [yf, Pc, Vm] = coreshell(p, q)
    v = p.coreRdivshR;
    sig = p.sig_fuzzy;
    Rg = p.Rg;
    if p.sigCore ==0
        r1 = v;
        nr1 = 1;
    else
        r1 = linspace(0.01, 0.99, 15);
        nr1 = schultz(v, p.sigCore*v, r1);
    end
    
    if p.sigshellR==0
        r2 = p.shellR;
        nr2 = 1;
    else
        r2 = linspace(1, p.shellR+3*p.sigshellR*p.shellR, 30);
        nr2 = schultz(p.shellR, p.sigshellR*p.shellR, r2);
    end

    Vchain = 4*pi/3*(Rg).^3;
    Pc = p.Nc*(p.delta_rho_chain*Vchain)^2*Pchain(q, Rg);
    if v == 0
        eden = [p.sh_rho, p.solvent_rho];
    else
        eden = [p.core_rho, p.sh_rho, p.solvent_rho];
    end
    Vm = 0;
    yf = zeros(size(q));
    for mm=1:numel(nr2)
        for kk=1:numel(nr1)
            if v == 0
                Ro = r2(mm);
            else
                Rc = r1(kk)*r2(mm);
                Ro = r2(mm);
            end
            if v ~= 0
                Vo = 4*pi/3*Ro^3;
                Vc = 4*pi/3*Rc.^3;
                P = ((eden(2)-eden(3))*Vo*Ps(q, Ro).*exp(-(q*sig).^2/2) ...
                    + (eden(1)-eden(2))*Vc*Ps(q, Rc)).^2;
            else
                Vo = 4*pi/3*Ro^3;
                P = ((eden(end)-eden(end-1))*Vo*Ps(q, Ro).*exp(-(q*sig).^2/2)).^2;
            end
            Vm = Vm + Vo*nr2(mm)*nr1(kk);
            yf = yf + P(:)*nr2(mm)*nr1(kk);
        end
    end
    
function y = Ps(q, R)
    y = 3*(sin(R.*q) - R.*q.*(cos(R.*q)))./((R.*q+eps).^3);
function y = Pchain(q, Rg)
    x = (q*Rg).^2;
    y = 2*(exp(-x)-1+x)./x.^2;
    
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