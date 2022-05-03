function [out, report] = FitLee_core_fuzzycorona_sphere(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Core-corona micelle model, where corona is fuzzy. ' ,...
'I(q) = fn0*(Sq*P(q; r0, sig0)) + background',...
'    background = poly1*q.^poly2 + poly3*q + poly4 + SF_userBG*UsersBackground',...
'    P(q) = fn0([(rho_shell-rho_solv)*Vo*Ps(qRo)exp(-(q*sig)^2/2) -',...
'           - (rho_shell-rho_core)*Vc*Ps(qRo*v)]^2',...
'           + Nc*(delta_rhochain*Vch)^2*Pc(q*Rg))',...
'Parameters',...
'  fn0 : volume fraction of particle (A^-3)',...
'  rho_shell : electron density of shell (A^-3)',...
'  rho_core : electron density of core (A^-3)',...
'  delta_rhochain : Excess electron density of chain (A^-3)',...
'  Ro : Outer radius (shell) (A)',...
'  v : ratio of radii of core and shell (v = Rc/Rs)',...
'  sig: Density gradient at shell/solvent interface (A)',...
'  Nc : Number chains per particle',...
'  Rg : Rg of the chain (A)',...
'  D : Interparticle distance (A) (hard sphere potential S(q))',...
'  vf : volume fraction of particle 0 (hard sphere potential S(q)) ',...
'  SF_userBG : Scale factor for the user input background',...
'  UsersBackground : Two column formatted user"" input background',...
'  ',...
'',...
'Note that',...
'I(q)_diff_cross_section = (delta_rho*r_e)^2*f_n*P(q)',...
'         , where P(0) = V_p^2'...
' ref: http://aip.scitation.org/doi/pdf/10.1063/1.1493771',...
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