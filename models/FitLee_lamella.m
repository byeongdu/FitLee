function [out, report] = FitLee_lamella(varargin)
% Usage:
% FitLee_lamella(N of peaks)
% FitLee('FitLee_multivoigt', bestP)
%
% Byeongdu Lee

FitLee_helpstr = {'SAXS from multilayer lamella (default 4 layers)' ,...
'I(q) = I0*I(q, thick, density, domainsize, microstrain, DW)/q^n + background',...
'    I(q) = lamella(q, thick, density, domainsize, microstrain, DW)',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'Parameters',...
'  I0 : Scaling factor',...
'  thick : An array of lamella thicknesses (A) ',...
'  density : An array of electron densities of the lamella (/A^3)',...
'  domainsize: Domain size of lamella stack, determing peak width (A).',...
'  microstrain : Micro strain of the structure, determing peak width (frac.)',...
'  DW : Debye-Waller factor, decaying intensities as q goes high (frac.)',...
'  ispowder : Is this a powder sample? 1 then n = 2 otherwise n = 0',...
'  ',...
'Byeongdu Lee (blee@anl.gov)'};

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
    p = 4;
    isini = 1;
end

%% initialize fit parameter bestP
if isini
    Nf = p;
    bestP = [];

    bestP.I0 = 1;
    for i=1:Nf
        bestP.(['thick', num2str(i)]) = 100;
        bestP.(['e_density', num2str(i)]) = rand(1);
    end
    bestP.DW = 0.005;
    bestP.domainsize = 2000;
    bestP.microstrain = 0.02;
    bestP.ispowder = 1;
    % Need 4 parameters for background.
    bestP.poly1 = 1E-8;
    bestP.poly2 = -2;
    bestP.poly3 = 0;
    bestP.poly4 = 0.001;
    assignin('base', 'bestP', bestP);
    out = bestP;
    return
end

%% fitting codes .......................
%        back = zeros(size(x));
%        if isfield(p, 'poly')
%            back = p.poly(1)*q.^p.poly(2);
%            for i=numel(p.poly):-1:3
%                back = back + p.poly(i)*q.^(numel(p.poly)-i);
%            end
%        end
fn = fieldnames(p);
Nf = (numel(fn) - 9)/2;
thick = zeros(size(Nf));
density = zeros(size(Nf));
for i=1:Nf
    thick(i) = p.(sprintf('%s%i', 'thick', i));
    density(i) = p.(sprintf('%s%i', 'e_density', i));
end
% lamella(q, thick, density, domainsize, microstrain, DW, ispowder)
[Iq, Fq, Fq2] = lamella(q, thick, density, p.domainsize, p.microstrain, p.DW, p.ispowder);
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
out = p.I0*Iq + back;

if nargout == 2
    %x = 0:1:(max(p.r0, p.r1)+max(p.sig0, p.sig1)*10);
    fprintf('MillerIndex   F(qhkl)    I(qhkl, DW)\n');    
    for k=1:Nf
        fprintf('h        %0.3e        #0.3e\n', k, Fq(k), Fq2(k));
    end    
    report = '';
end

%varargout{2} = fit;
%err = sum( (y-predY).^2);

