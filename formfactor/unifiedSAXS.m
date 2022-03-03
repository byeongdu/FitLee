function [out, report] = unifiedSAXS(varargin)
global FitLee_helpstr
FitLee_helpstr = {'Unified SAXS by Greg Beaucage. ' ,...
' P(q) = G*exp(-x.^2*Rg^2/3) + B*(erf(x*Rg/sqrt(6)).^3./x).^P',...
'   I0 : I0',...
'   B : Porod intensity',...
'   Rg : Rg', ...
'   P : Porod exponent',...
'I(q) = P(q)*Sq + background',...
'    Sq = powI*q.^PorodExp + S(q; D, vf)',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'Parameters',...
'  D : Interparticle distance (A) (hard sphere potential S(q))',...
'  vf : volume fraction of particle 0 (hard sphere potential S(q)) ',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    1. G. Beaucage, J. appl. Cryst(1995), 28, 717-728. '};

if numel(varargin) > 1
    p = varargin{1};
    q = varargin{2};
    isini = 0;
elseif numel(varargin) == 1
    p = varargin{1};
    isini = 1;
elseif numel(varargin) == 0
    p = 1;
    isini = 1;
end

%% initialize fit parameter bestP
if isini
    Nf = p;
    bestP = [];
    bestP.I0 = 10000;
    bestP.B = 1E-6;
    bestP.Rg = 300;
    bestP.P = 4;
    bestP.D = 100;
    bestP.vf = 0.0;
    bestP.powI = 1E-4;
    bestP.PorodExp = -4;
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
        error('This function is for fitting a set of data, for now')
    end
    q = q{1};
end

q = q(:);

%Pq1 = saxs_poly_fractalparticle(q, p.I0, p.r0, p.sig0, p.Df);
Pq = p.I0*(exp(-q.^2*p.Rg^2/3) + p.B*(erf(q*p.Rg/sqrt(6)).^3./q).^p.P);

Sq1 = strfactor2(q, p.D, p.vf);
Sq = p.powI*q.^p.PorodExp + Sq1;
%Sq = p(4)*q(:).^p(5) + strfactor_2Dpara(q, p(6), p(7));
%Iq = p(1)*Imat*nr1/sum(nr1)
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
%pnumberfraction = p.f0;
out = Pq.*Sq + +back;
if isnan(out)
    out = ones(size(out));
end



% function [y, name, pnames, pin] = unifiedSAXS_original(x, p, flag)
% % [y, name, pnames, pin] = unifiedSAXS(x, p, flag)
% % Unified model : G. Beaucage, J. appl. Cryst(1995), 28, 717-728.
% % y = G*exp(-x.^2*Rg^2/3) + B*(erf(x*Rg/sqrt(6)).^3./x).^P;
% %   G  = p(1);
% %   B = p(2);
% %   Rg = p(3);
% %   P = p(4);
%  
% if nargin == 2;
%     if isstruct(p)
%         G = p.G;
%         B = p.B;
%         Rg = p.Rg;
%         P = p.P;
%     else
%         G  = p(1);
%         B = p(2);
%         Rg = p(3);
%         P = p(4);
%     end
% %y=1
% else
% 	y=[];
% 	name='Unified SAXS';
% 	pnames=str2mat('G','B', 'Rg', 'P');
% 	if flag==1, pin=[100 0.0005 500 2]; else pin = p; end
% end
% end