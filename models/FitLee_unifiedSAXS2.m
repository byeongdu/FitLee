function [out, report] = FitLee_unifiedSAXS2(varargin)
FitLee_helpstr = {'Unified SAXS by Greg Beaucage. ' ,...
' P(q) = G*exp(-x.^2*Rg^2/3) + B*(erf(x*Rg/sqrt(6)).^3./x).^P',...
'   I0 : I0',...
'   B : Porod intensity',...
'   Rg : Rg', ...
'   P : Porod exponent',...
' P1(q) = guinierporodmodel(q, G, d, Rg, s)',...
'   I0 : Guinier I0',...
'   Rg : Rg of a particle',...
'   P : Porod exponent',...
'   s : dimension of the particle',... 
'I(q) = powI*q.^PorodExp + P1(q) + P0(q)*Sq + background',...
'    Sq = S(q; D, vf)',...
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
    bestP.I0_1 = 10000;
    bestP.s_1 = 1E-6;
    bestP.Rg_1 = 300;
    bestP.P_1 = 4;
    bestP.I0_0 = 10000;
    bestP.B_0 = 1E-6;
    bestP.Rg_0 = 300;
    bestP.P_0 = 4;
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

q = abs(q(:));

%Pq1 = saxs_poly_fractalparticle(q, p.I0, p.r0, p.sig0, p.Df);
Pq_0 = p.I0_0*(exp(-q.^2*p.Rg_0^2/3) + p.B_0*(erf(q*p.Rg_0/sqrt(6)).^3./q).^p.P_0);
Pq_1 = guinierporodmodel(q, p.I0_1, p.P_1, p.Rg_1, p.s_1);
Sq1 = strfactor2(q, p.D, p.vf);
%Sq = p.powI*q.^p.PorodExp + Sq1;
%Sq = p(4)*q(:).^p(5) + strfactor_2Dpara(q, p(6), p(7));
%Iq = p(1)*Imat*nr1/sum(nr1)
back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
%pnumberfraction = p.f0;
out = p.powI*q.^p.PorodExp + Pq_0.*Sq1 + Pq_1 + back;
if isnan(out)
    out = ones(size(out));
end

if nargout == 2
    %q = linspace(1E-5, 1, 2);
    Q0 = invariant(q, Pq_0, [q(1), q(end)]);
    Q1 = invariant(q, Pq_1, [q(1), q(end)]);
    fprintf('Invariant for the unified fit : %0.3e cm^{-1}A^{-3}.\n', Q0);
    fprintf('Invariant for the Guiner-Porod fit : %0.3e cm^{-1}A^{-3}.\n', Q1);
    report = '';
end

