function [out, report] = FitLee_cluster_coreshellfuzzysphere(varargin)
FitLee_helpstr = {'Polydisperse mass-fractal particles. ' ,...
'I(q) = I0*P(q; avgR, z, Df)*Sq + background',...
'    Sq = guinierporodmodel(q, I0_2, s_2, Rg_2, P_2) + S(q; D, vf)',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'Parameters',...
' P(q)',...
'  I0 : Scale factor I0',...
'  r0 : Schultz size distribution, radius peak (A)',...
'  sig0: Schultz size distribution, FWHM (A)',...
'  Df : mass fractal dimension.',...
' S(q) = P1(q) + S_inf(q)',...
' P1(q) = guinierporodmodel(q, I0, d, Rg, s)',...
'   I0 : Guinier I0, meaning the number of primary particles forming a cluster',...
'   Rg : Rg of the cluster',...
'   P : Porod exponent of the cluster',...
'   s : dimension of the cluster',... 
' S_inf(q)',...
'  D : Interparticle distance (A) (hard sphere potential S(q))',...
'  vf : volume fraction of particle 0 (hard sphere potential S(q)) ',...
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
% primary particle    
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


    % structure factor
    bestP.I0_Cluster = 1;
    bestP.s_Cluster = 2;
    bestP.Rg_Cluster = 300;
    bestP.P_Cluster = 4;
    bestP.D_Cluster = 100;
    bestP.vf_Cluster = 0.0;
    bestP.powI_Cluster = 1E-8;
    bestP.PorodExp_Cluster = -4;
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
if iscell(q)
    if numel(q) > 1
        error('FitLee_schultzsphere.m is for fitting a set of data, for now')
    end
    q = q{1};
end

q = q(:);

    p2.fn0 = 1;
    p2.core_rho = p.core_rho;
    p2.coreRdivshR = p.coreRdivshR;
    p2.sigCore = p.sigCore;
    
    p2.sh_rho = p.sh_rho;
    p2.shellR = p.shellR;
    p2.sigshellR = p.sigshellR;
    p2.sig_fuzzy = p.sig_fuzzy;
    
    p2.solvent_rho = p.solvent_rho;
    
    p2.Nc = p.Nc;
    p2.delta_rho_chain = p.delta_rho_chain;
    p2.Rg = p.Rg;
    
    p2.D = p.D;
    p2.vf = p.vf;

    p2.poly1 = 0;
    p2.poly2 = 0;
    p2.poly3 = 0;
    p2.poly4 = 0;
    p2.SF_userBG = 0;
    pqr = FitLee_core_fuzzycorona_sphere(p2, q);
    Pq1 = pqr(:,1);
    Pq1 = Pq1*p.I0_PP;

Sq0 = strfactor2(q, p.D_Cluster, p.vf_Cluster);
Sq1 = strfactor2(q, p.D, p.vf);
Pq_1 = guinierporodmodel(q, p.I0_Cluster, p.P_Cluster, p.Rg_Cluster, p.s_Cluster);
Poq = p.powI_Cluster*q.^p.PorodExp_Cluster;

Sq = Poq + Pq_1.*Sq0 + Sq1;
Iq = Pq1(:).*Sq(:);

%Sq = p(4)*q(:).^p(5) + strfactor_2Dpara(q, p(6), p(7));
%Iq = p(1)*Imat*nr1/sum(nr1)
back = p.poly1*q(:).^p.poly2 + p.poly3*q(:) + p.poly4;
%pnumberfraction = p.f0;
Iqa = Iq + back;
out = [Iqa(:), Iq(:), Pq1(:), back(:)];
if isnan(out)
    out = ones(size(out));
end

