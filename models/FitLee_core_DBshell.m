function [out, report] = FitLee_core_DBshell(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Core-shell model, where the shell is a DB type. ' ,...
'I(q) = fn0*(Sq*P(q; r0, sig0)) + background',...
'    background = poly1*q.^poly2 + poly3*q + poly4 + SF_userBG*UsersBackground',...
'    P(q) = |(rhoc-rhosh)Fsp(Rc)+(rhosh-rhosolv)F_DBshell|^2',...
'         ~ |(rhoc-rhosh)|^2|Fsp(Rc)|^2+|(rhosh-rhosolv)F_DBshell|^2',...
'Parameters',...
'  fn0 : volume fraction of particle ',...
'  rho_shell : electron density of shell (A^-3)',...
'  rho_core : electron density of core (A^-3)',...
'  xsi : Correlation length of the shell (A)',...
'  sig: Density gradient at shell/solvent interface (A)',...
'  D : Interparticle distance (A) (hard sphere potential S(q))',...
'  vf : volume fraction of particle 0 (hard sphere potential S(q)) ',...
'  ',...
'',...
'Note that',...
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
    bestP.core_rho = 4.66;
    bestP.Rcore = 100;
    bestP.sigcore = 10;
    
    bestP.sh_rho = 0.34;
    bestP.sh_Rg = 150;
    bestP.chainRg = 50;
    bestP.Nchains = 1000;
    bestP.polymerMw = 2000;
    bestP.monomerMw = 104.15;
    bestP.polymerdensity = 1.04;
    
    bestP.solvent_rho = 0.284;
    
    bestP.D = 200;
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

[y, V2, V, ~, ~, Fq] = SchultzsphereFun(q, p.Rcore, p.sigcore);
Isp = y*V2;
pnumberfraction = p.fn0/(V);
f = pnumberfraction*r_e^2;
Isp = f*(p.core_rho - p.sh_rho)^2*Isp;

delta_rho =  p.sh_rho - p.solvent_rho;
xsi = p.sh_Rg/sqrt(6); % see my chem rev. text below eqn 63.
VDB = 8*pi*xsi^3; % V = invariant of DB.
volf = VDB;
%volf = f*(p.sh_Rg*sqrt(5/3))^3*4*pi/3*5.2705;
yDB = DebyeBueche(q, [delta_rho, volf, xsi, 0]);

%Icross = f*(p.core_rho - p.sh_rho)*V*exp(-q.^2*3/5*p.Rcore^2/3/10).*sqrt(yDB);
yDB = f*yDB;

% chain
monomervol = p.monomerMw/p.polymerdensity/(6.02E23)/Angstrom2Centimeter^3; % unit cm^3
Nmonomer = p.polymerMw/p.monomerMw;
yP = delta_rho^2*monomervol^2*Nmonomer^2*Pchain(q, p.chainRg); % scattering from a single chain.
yP = f*p.Nchains*yP;

Sq = strfactor2(q, p.D, p.vf);

back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
yDB = yDB/Angstrom2Centimeter;
yP = yP/Angstrom2Centimeter;
Isp = Isp/Angstrom2Centimeter;
Iq0 = (yDB + yP + Isp).*Sq;
out = Iq0 + back;
out = [out(:), Isp(:), yDB(:), yP(:), back(:)];
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
function y = Pchain(q, Rg)
    x = (q*Rg).^2;
    y = 2*(exp(-x)-1+x)./x.^2;

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