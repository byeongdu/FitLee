function out = FitLee_lamellaparacrystal(varargin)
    % Paracrystal theory equation for triblock copolymer.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'diblock copolymer diffraction analysis. ' ,...
'Based on Paracrystal theory',...
'',...
'Lamella structure',...
'deltarho : electron density difference between C and A phases',...
'I0 : Number of lamella',...
'dc : thickness(A) of C lamella',...
'da : thickness(A) of A lamella',...
'sigmadc : std thickness(A) of C lamella',...
'sigmada : std thickness(A) of A lamella',...
'',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    1. R.J. Roe book '};

if numel(varargin) > 1
    p = varargin{1};
    q = varargin{2};
    isini = 0;
elseif numel(varargin) == 1
    p = varargin{1};
    isini = 1;
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
    bestP.I0 = 1;
    bestP.deltarho = 2.0;
    bestP.da = 110;
    bestP.sigmada = 40;
    bestP.dc = 95;
    bestP.sigmadc = 20;
    % Need 4 parameters for background.
    bestP.poly1 = 0;
    bestP.poly2 = -2;
    bestP.poly3 = 0;
    bestP.poly4 = 0;
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
q = q(:);

N = p.I0;
deltarho = p.deltarho;
da = p.da;
dc = p.dc;
sigmaa = p.sigmada;
sigmac = p.sigmadc;

Iq = lamella2(q,[N, deltarho, da, dc, sigmaa, sigmac]);
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
out = Iq + back;

if isnan(out)
    out = ones(size(out))*1E20;
end

end
