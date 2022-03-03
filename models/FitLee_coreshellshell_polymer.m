function [out, report] = FitLee_coreshellshell_polymer(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Schultz polydisperse sphere fit. ' ,...
'I(q) = I0*(Sq*P(q; r0, sig0) + I1*P(q; r1, sig1)) + background',...
'    Sq = powI*q.^PorodExp + S(q; D, vf)',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'CF. if you have data in absolute unit, consider use FitLee_schultzsphere2',...
'',...
'Byeongdu Lee (blee@anl.gov)',...
'Ref: ',...
'    1. Kwon et al. Nature Materials, 2015, 14, 215–223. ',...
'    2. Wang et al. J. Phys. Chem. C. 2013. 117(44), 22627. '};

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
    bestP.D0 = 100;
    bestP.vf0 = 0.01;
    bestP.I0 = 4.6095e-13;
    bestP.coreRdivshR = 0.245672;
    bestP.sigCore = 0.097654;
    bestP.shellR = 247.691904;
    bestP.sigshellR = 0.012843;
    bestP.sh2thick = 42.636334;
    bestP.sigsh2thick = 0.046388;
    bestP.core_rho = 0.38;
    bestP.sh1_rho = 0.39;
    bestP.sh2_rho = 0.42;
    bestP.solvent_rho = 0.3344;
    bestP.d_rho_chain = 1;
    bestP.Nc = 100;
    bestP.Rg = 50;
    % Need 4 parameters for background.
    bestP.I1 = 0;
    bestP.r1 = 35;
    bestP.sig1 = 3.5;
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
if iscell(q)
    if numel(q) > 1
        error('FitLee_schultzsphere.m is for fitting a set of data, for now')
    end
    q = q{1};
end

q = q(:);

method = 1;

Rg = p.Rg;
Vchain = 4*pi/3*(Rg).^3;
Pc = p.Nc*(p.d_rho_chain*Vchain)^2*Pchain(q, Rg);


Iq = coreshell(p, q, method);

[Pq1, V1] = SchultzsphereFun(q, p.r1, p.sig1);
Pq1 = V1*Pq1(:);
Sq0 = strfactor2(q, p.D0, p.vf0);
Sq1 = strfactor2(q, p.D, p.vf);
%Sq = p.powI*q.^p.PorodExp + Sq1;

back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;

Iq1 = p.I0*Iq.*Sq0;
Iq2 = p.I1*Pq1.*Sq1;
Pc = Pc.*Sq1;
out = Iq1 + Iq2 + Pc + back;
out = [out(:), Iq1(:), Iq2(:), Pc(:), back(:)];

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

function y = Pchain(q, Rg)
    x = (q*Rg).^2;
    y = 2*(exp(-x)-1+x)./x.^2;

function Iq = coreshell(p, q, method)

% method = 1;

if method == 1
    if p.sigCore ==0
        r1 = p.coreRdivshR;
        nr1 = 1;
    else
        r1 = linspace(0.01, 0.99, 15);
        %nr1 = schultz(p.coreRdivshR, p.sigCore, r1);
        nr1 = schultz(p.coreRdivshR, p.sigCore*p.coreRdivshR, r1);
    end
    if p.sigshellR==0
        r2 = p.shellR;
        nr2 = 1;
    else
        %r2 = linspace(1, p.shellR+3*p.sigshellR, 30);
        r2 = linspace(1, p.shellR+3*p.sigshellR*p.shellR, 30);
        %nr2 = schultz(p.shellR, p.sigshellR, r2);
        nr2 = schultz(p.shellR, p.sigshellR*p.shellR, r2);
    %Rmat.r2 = r2;
    end
%     Rmat.r2 = r2;
%     Rmat.dist2 = nr2;
    if p.sigsh2thick ==0
        r3 = p.sh2thick;
        nr3 = 1;
    else
%         r3 = linspace(0.1, p.sh2thick+3*p.sigsh2thick, 30);
%         nr3 = schultz(p.sh2thick, p.sigsh2thick, r3);
        r3 = linspace(0.1, p.sh2thick+3*p.sigsh2thick*p.sh2thick, 30);
        nr3 = schultz(p.sh2thick, p.sigsh2thick*p.sh2thick, r3);
    end
%     Rmat.r3 = r3;
%     Rmat.dist3 = nr3;

    yf = zeros(size(q));
    for mm=1:numel(nr2);
        for kk=1:numel(nr1)
            rad = [];
            for ll=1:numel(nr3)
                if p.coreRdivshR == 0
                    rad = [rad;r2(mm),r2(mm) + r3(ll)];
                else
                    rad = [rad;r1(kk)*r2(mm), r2(mm), r2(mm)+r3(ll)];
                end
            end
            if p.coreRdivshR == 0
                eden = [p.sh1_rho, p.sh2_rho, p.solvent_rho];
            else
                eden = [p.core_rho, p.sh1_rho, p.sh2_rho, p.solvent_rho];
            end
            yt = multilayersphere2(q, rad,eden, nr3);
            yf = yf + yt(:)*nr2(mm)*nr1(kk);
        end
    end
end
if method == 2
    Rmat = evalin('base', 'Rmat');

    if p.coreRdivshR ==0
        nr1 = 1;
        nr2 = schultz(p.shellR, p.sigshellR, r2);
        %nr2 = 1;
        %Rmat.dist2 = 1;
        %Rmat.dist2 = nr2;
        if p.sigsh2thick ==0
            nr3 = 1;
        else
            nr3 = schultz(p.shellR+p.sh2thick, p.sigsh2thick, r2);
        end
        %dist3 = nr3;
        eden = [p.sh1_rho, p.sh2_rho, p.solvent_rho];
        dis = {nr2, nr3};
        Fm = {Rmat.Fmat2, Rmat.Fmat3};
    else
        eden = [p.core_rho, p.sh1_rho, p.sh2_rho, p.solvent_rho];
        nr1 = schultz(p.coreRdivshR*p.shellR, p.coreRdivshR*p.sigCore*p.shellR, Rmat.r1);
        Rmat.dist1 = nr1;
        nr2 = schultz(p.shellR, p.sigshellR, Rmat.r2);
        Rmat.dist2 = nr2;
        nr3 = schultz(p.sh2thick+p.shellR, p.sigsh2thick, Rmat.r2);
        Rmat.dist3 = nr3;
        dis = {nr1, nr2, nr3};
        Fm = {Rmat.Fmat1, Rmat.Fmat2, Rmat.Fmat3};
    end
    yf = multisphere(Fm, eden, dis);
    assignin('base', 'Rmat', Rmat)
end

Iq = yf(:)/sum(nr1)/sum(nr2)/sum(nr3);

function Pq = multisphere(Fm, eden, dis)
eden = eden(1:end-1) - eden(2:end);
Fmat1 = Fm{1};dist1 = dis{1};
Fmat2 = Fm{2};dist2 = dis{2};
if numel(Fm) == 3
    Fmat3 = Fm{3};
    dist3 = dis{3};
end
    
numq = size(Fmat1, 1);
if numel(Fm) == 2
    for i=1:numq
        [f1,f2]=meshgrid(Fmat1(i,:),Fmat2(i,:));
        [d1,d2] = meshgrid(dist1, dist2);
        f = f1*eden(1)+f2*eden(2);

        f = abs(f).^2.*d1.*d2;
        f = sum(sum(f));
        Pq(i) = f;
    end
end
if numel(Fm) == 3
    for i=1:numq
        [f1,f2,f3]=meshgrid(Fmat1(i,:),Fmat2(i,:),Fmat3(i,:));
        [d1,d2,d3] = meshgrid(dist1, dist2, dist3);
        f = f1*eden(1)+f2*eden(2)+f3*eden(3);
        f = abs(f).^2.*d1.*d2.*d3;
        f = sum(sum(sum(f)));
        Pq(i) = f;
    end
end