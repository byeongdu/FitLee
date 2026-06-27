function [out, report] = FitLee_core_fuzzyshell(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Core + fuzzy shell (erf-tapered), absolute unit $(cm^{-1}). $' ,...
'',...
'$I(q) = (fn0/\left<V_m\right>) \cdot r_e^2 \cdot P_{cs}(q) \cdot S(q; D, v_f) / 10^{-8} + I_b$',...
'$\qquad    P_{cs}(q) = \left| (\rho_{core} - \rho_{solv}) V_c P_s(q R_c) + (\rho_{sh} - \rho_{solv}) V_{shell\_eff}(q) \right|^2$',...
'$\qquad    $ The shell density profile decays as $\frac{1}{2}\,\mathrm{erfc}((r - R_c - t_{sh}) / (\sqrt{2}\sigma_{fuzzy}))$',...
'$\qquad    \left<V_m\right>$ from schultzRg$(R_{core}, sig_{core})$ (mean particle volume)',...
'$\qquad    I_b = poly1 \cdot q^{poly2} + poly3 \cdot q + poly4$',...
'',...
'Note: no chain term in this model (vs FitLee_core_fuzzycorona_sphere which adds one).',...
'The shell is described by an erf-tapered density profile of thickness $t_{sh}$ with',...
'Gaussian roughness $\sigma_{fuzzy}$ at the outer interface.',...
'',...
'$Parameters$',...
'$\quad  fn0$ : volume fraction of particles (dimensionless)',...
'$\quad  core\_rho, sh\_rho, solvent\_rho$ : electron densities $(A^{-3})$',...
'$\quad  R_{core}, sig_{core}$ : Schultz radius peak (A) and FWHM (A) for the core',...
'$\quad  sh_{thick}$ : nominal shell thickness (A)',...
'$\quad  sig\_fuzzy$ : Gaussian roughness $\sigma$ (A) at shell/solvent interface',...
'$\quad  D, vf$ : hard-sphere $S(q)$ parameters',...
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
    bestP.Rcore = 0.5;
    bestP.sigcore = 0.1;
    
    bestP.sh_rho = 0.4;
    bestP.sh_thick = 50;
    bestP.sig_fuzzy = 10;
    
    bestP.solvent_rho = 0.3344;
    
    bestP.D = 100;
    bestP.vf = 0.1;

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



r_e = 2.818E-5; % Angstrom
Angstrom2Centimeter = 1E-8;
%r_e2 = (r_e*Angstrom2Centimeter)^2;

q = q(:);
P = coreshell(p, q);
Sq = strfactor2(q, p.D, p.vf);

back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
[~, Vm] = schultzRg(p.Rcore, p.sigcore);
pnumberfraction = p.fn0/Vm/Angstrom2Centimeter;
f = pnumberfraction*r_e^2;
Iq0 = f*P.*Sq;
out = Iq0 + back;
out = [out(:), Iq0(:), back(:)];
if isnan(out)
    out = ones(size(out));
end
if nargout == 2
    x = 0:1:(p.Rcore + p.sigcore*10);
    nr0 = schultz(p.Rcore, p.sigcore, x);
    nr = nr0;
    figure;
    subplot(1,2,1)
    plot(x, nr);xlabel('Radius (A)');ylabel('n(r)')
    subplot(1,2,2)
    r = 0:1:(p.Rcore+p.sh_thick+p.sig_fuzzy*5/2);
    rho = (p.sh_rho - p.solvent_rho)*(1/2-1/2*erf((r-p.Rcore-p.sh_thick)/(sqrt(2)*p.sig_fuzzy)));
    rho(r<p.Rcore) = p.core_rho-p.solvent_rho;
    plot(r, rho)
    report = '';
end
function yf = coreshell(p, q)
    sig = p.sig_fuzzy;
    if p.sigcore ==0
        r1 = p.Rcore;
        nr1 = 1;
    else
        [nr1, r1] = schultzdist99(p.Rcore, p.sigcore, 15);
    end
    
    eden = [p.core_rho, p.sh_rho, p.solvent_rho];
    yf = zeros(size(q));
    for kk=1:numel(nr1)
        Rc = r1(kk);
        Rsh = Rc + p.sh_thick;
        P = (eden(2)-eden(3))*sphereamp(q, Rsh).*exp(-(q*sig).^2) ...
            + (eden(1)-eden(2))*sphereamp(q, Rc);
        yf = yf + P(:).^2*nr1(kk);
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