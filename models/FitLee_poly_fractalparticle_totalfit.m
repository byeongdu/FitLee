function [out, report] = FitLee_poly_fractalparticle_totalfit(varargin)
FitLee_helpstr = {'Polydisperse mass-fractal particles, Porod tail, with up to 4 external backgrounds.' ,...
'',...
'$I(q) = I_{0,PP} \cdot P_{PP}(q) \cdot S(q) + I_{pw} + I_b$',...
'$\qquad    P_{PP}(q) = $ saxs_poly_fractalparticle$(q, 1, r_{0,PP}, sig_{0,PP}, Df_{PP})$ if $Df_{PP} < 3$',...
'$\qquad    P_{PP}(q) = $ SchultzsphereFun$(q, r_{0,PP}, sig_{0,PP})$ otherwise',...
'$\qquad    I_c(q) = P_{cluster}(q) \cdot S_{HS}(q; R_{h,Cluster}, v_{f,Cluster})$',...
'$\qquad    S_1(q) = S_{HS}(q; R_{h,PP}, v_{f,PP})$',...
'$\qquad    S(q) = \max(I_c, S_1) + \dfrac{1}{(|I_c - S_1| + 1.5)^4} - \dfrac{1}{(2.5)^4}$',...
'$\qquad    P_{cluster}(q) = $ guinierporodmodel$(q, I_{0,Cluster}, P_{Cluster}, Rg_{Cluster}, s_{Cluster})$',...
'$\qquad    I_{pw} = powI_{Cluster} \cdot q^{PorodExp_{Cluster}}$',...
'$\qquad    I_b = poly1 \cdot back_1 + poly2 \cdot back_2 + poly3 \cdot back_3 + poly4 \cdot back_4$',...
'',...
'Note: $poly1..poly4$ are weights on UP TO FOUR user-loaded backgrounds named',...
'$back1, back2, back3, back4$ in the base workspace (each a two-column $q$-vs-$I$ array,',...
'interpolated to the fit $q$). Missing backgrounds default to 0. This is NOT the',...
'polynomial form used elsewhere.',...
'',...
'$Parameters$',...
'$\quad  I0\_PP$ : scale factor for the primary fractal particle',...
'$\quad  r_{0,PP}, sig_{0,PP}$ : Schultz radius peak (A) and FWHM (A)',...
'$\quad  Df_{PP}$ : mass-fractal dimension (use $< 3$ for fractal, $\geq 3$ for Schultz sphere)',...
'$\quad  Rh\_PP, vf\_PP$ : hard-sphere $S\_1(q)$ for PP arrangement',...
'$\quad  I0\_Cluster, P\_Cluster, Rg\_Cluster, s\_Cluster$ : Guinier-Porod for cluster',...
'$\quad  Rh\_Cluster, vf\_Cluster$ : hard-sphere $S(q)$ folded into $I\_c(q)$',...
'$\quad  powI\_Cluster, PorodExp\_Cluster$ : additive low-q power-law $I_{pw}$',...
'$\quad  poly1..poly4$ : weights on user backgrounds $back1..back4$ from base workspace',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    1. R. Besselink and J. E. ten Elshof, J. Appl. Cryst., (2014), 47, 1606. '};

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
    bestP.I0_PP = 1;
    bestP.r0_PP = 100;
    bestP.sig0_PP = 10;
    bestP.Df_PP = 2;
    bestP.Rh_PP = 100;
    bestP.vf_PP = 0.0;
    % structure factor
    bestP.I0_Cluster = 1;
    bestP.s_Cluster = 2;
    bestP.Rg_Cluster = 300;
    bestP.P_Cluster = 4;
    bestP.Rh_Cluster = 100;
    bestP.vf_Cluster = 0.0;
    bestP.powI_Cluster = 1E-8;
    bestP.PorodExp_Cluster = -4;
    % Need 4 parameters for background.
    bestP.poly1 = 1;
    bestP.poly2 = 1;
    bestP.poly3 = 1;
    bestP.poly4 = 1;
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
    back = evalin('base', 'back1');
    back1 = interp1(back(:,1), back(:,2), q);
catch
    back1 = 0;
end
try
    back = evalin('base', 'back2');
    back2 = interp1(back(:,1), back(:,2), q);
catch
    back2 = 0;
end
try
    back = evalin('base', 'back3');
    back3 = interp1(back(:,1), back(:,2), q);
catch
    back3 = 0;
end
try
    back = evalin('base', 'back4');
    back4 = interp1(back(:,1), back(:,2), q);
catch
    back4 = 0;
end

q = q(:);
RgPP = [];
if p.Df_PP < 3
    [Pq1, RgPP] = saxs_poly_fractalparticle(q, 1, p.r0_PP, p.sig0_PP, p.Df_PP);
else
    Pq1 = SchultzsphereFun(q, p.r0_PP, p.sig0_PP);
    Pq1 = Pq1*p.I0_PP;
%    Pq0 = V1*Pq1(:);
end

Sq0 = strfactor2(q, p.Rh_Cluster, p.vf_Cluster);
Sq1 = strfactor2(q, p.Rh_PP, p.vf_PP);
Pq_1 = guinierporodmodel(q, p.I0_Cluster, p.P_Cluster, p.Rg_Cluster, p.s_Cluster);
Poq = p.powI_Cluster*q.^p.PorodExp_Cluster;

Ic = Pq_1.*Sq0;
%if p.vf_PP > 0
    k2 = Ic < Sq1;
    ind = find(k2==1);
    if numel(ind)>1
    if q(ind(1))>(pi/p.Rh_PP)
        p.I0_PP = p.I0_PP*1E10;
    end
    end
    Sq = (Ic+Sq1);
    Sq(k2) = Sq1(k2);
    Sq(~k2) = Ic(~k2);
% else
%     Sq = Ic+Sq1;
% end
Sq = Sq + 1./(abs(Ic-Sq1)+1.5).^4-(1/(abs(1)+1.5)^4);

%k = q > 2*pi/p.D_PP;
%k = Ic2 < 1;
%Ic2(k) = 1;

%Sq = Poq + Sq;
Iq = p.I0_PP*Pq1(:).*Sq(:) + Poq;

%Sq = p(4)*q(:).^p(5) + strfactor_2Dpara(q, p(6), p(7));
%Iq = p(1)*Imat*nr1/sum(nr1)
%back = p.poly1*q(:).^p.poly2 + p.poly3*q(:) + p.poly4;
back = p.poly1*back1 + p.poly2*back2 + p.poly3*back3 + p.poly4*back4; 
%pnumberfraction = p.f0;
if length(back) ==1
    back = back*ones(size(Iq));
end
out = [Iq + back, p.I0_PP*Pq1(:), Sq(:), back];
if isnan(out)
    out = ones(size(out));
end

if nargout == 2
    x = 0:1:((p.r0_PP+p.sig0_PP)*10);
    nr0 = schultz(p.r0_PP, p.sig0_PP, x);
    nr = nr0;
    figure;subplot(2,1,1)
    plot(x, nr);xlabel('Radius (A)');ylabel('n(r)')
    subplot(2,1,2)
    Vr = nr(:).*x(:).^3;
    plot(x, Vr);xlabel('Radius (A)');ylabel('V(r)')
    volumedistribution = [x(:), Vr];
    assignin('base', 'volumedistribution', volumedistribution);
    %[zRg, V, ~, S] = schultzRg(p.r0, p.sig0);
    if ~isempty(RgPP)
        zRg = RgPP;
    else
        zRg = 0;
    end
    fprintf('Rg of particle0 : %0.3fA.\n', zRg);
    %[~, ~, ~, S] = schultzRg(p.r0, p.sig0);
    %S = S/4/100; % 4 pi R(A) ^2 --> pi*R(nm)^2
    report = '';
end
