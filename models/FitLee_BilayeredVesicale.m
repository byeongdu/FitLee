function [out, report] = FitLee_BilayeredVesicale(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'FitLee_BilayeredVesicale. ' ,...
'fsig : fractional sig',...
'I(q) = I0*(I_cs(q; r0, fsig0) + I1*P(q; r1, fsig1)*Sq) + background',...
'    Sq = powI*q.^PorodExp + S(q; D, vf)',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'CF. if you have data in absolute unit, consider use FitLee_schultzsphere2',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
};

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
    bestP.I0 = 4.6095e-13;
    bestP.Rc = 30;
    bestP.fsigCore = 0.097654;
    bestP.t_h = 3;
    bestP.t_h_fsig = 0.1;
    bestP.t_t = 0.5;
    bestP.head_rho = 0.39;
    bestP.tail_rho = 0.26;
    bestP.solvent_rho = 0.3344;
    bestP.I1 = 0;
    bestP.r1 = 35;
    bestP.fsig1 = 0.1;
    bestP.D = 100;
    bestP.vf = 0.01;
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

method = 1;

Iq = coreshell(p, q, method);

[Pq1, V1] = SchultzsphereFun(q, p.r1, p.r1*p.fsig1);
Pq1 = V1*Pq1(:);
Sq1 = strfactor2(q, p.D, p.vf);
%Sq = p.powI*q.^p.PorodExp + Sq1;

back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
out = p.I0*Iq + p.I1*Pq1.*Sq1 + back;
if isnan(out)
    out = ones(size(out))*1E20;
end

if nargout == 2
    x = 0:1:(p.shellR + p.sh2thick + 10);
    rho = zeros(size(x));
    t = x < p.coreRdivshR*p.shellR;
    
    rho(t) = p.core_rho;
    t = (x >= p.coreRdivshR*p.shellR) & (x < p.shellR);
    rho(t) = p.sh1_rho;
    t = (x >= p.shellR) & (x < p.shellR+p.sh2thick);
    rho(t) = p.sh2_rho;
    t = (x >= p.shellR+p.sh2thick);
    rho(t) = p.solvent_rho;
    
    figure;
    plot(x, rho);xlabel('Radius (A)');ylabel('\rho (R)')
    report = '';
end



function Iq = coreshell(p, q, method)

% method = 1;

if method == 1
    if p.fsigCore ==0
        r1 = p.Rc;
        nr1 = 1;
    else
        rmax = p.Rc + 6*p.Rc*p.fsigCore;
        Npnt = (1-0.01)/max(q)*rmax;
        if Npnt > 150
            Npnt = 150;
        end
        r1 = linspace(0.01, 1, Npnt);
        r1 = r1*rmax;
        %nr1 = schultz(p.coreRdivshR, p.sigCore, r1);
        nr1 = schultz(p.Rc, p.fsigCore*p.Rc, r1);
    end
    if p.t_h_fsig==0
        th = p.t_h;
        nr2 = 1;
    else
        %r2 = linspace(1, p.shellR+3*p.sigshellR, 30);
        th = linspace(0.01, 1, 15);
        th = th*(p.t_h + 6*p.t_h*p.t_h_fsig);
        %nr2 = schultz(p.shellR, p.sigshellR, r2);
        nr2 = schultz(p.t_h, p.t_h*p.t_h_fsig, th);
    end
    nr3 = 1;
    
    yf = zeros(size(q));
    for mm=1:numel(nr2);
        for kk=1:numel(nr1)
            Aq = (p.solvent_rho - p.head_rho)*sphereamp(q, r1(kk)) +...
                (p.head_rho-p.tail_rho)*sphereamp(q, r1(kk) + th(mm)) +...
                (p.tail_rho-p.head_rho)*sphereamp(q, r1(kk) + th(mm) + p.t_t) +...
                (p.head_rho-p.solvent_rho)*sphereamp(q, r1(kk) + 2*th(mm) + p.t_t);
            yf = yf + Aq(:).^2*nr2(mm)*nr1(kk);
        end
    end
end

Iq = yf(:)/sum(nr1)/sum(nr2)/sum(nr3);