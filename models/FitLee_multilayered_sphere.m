function [out, report] = FitLee_multilayered_sphere(varargin)
    % schultz polydisperse sphere with structure factors.
    % This function is only for fitting a set of data at a time.
    % if numel(varargin) == 1 and p is not a struct but a number, then
    % it generate set of default parameters for using this function.
FitLee_helpstr = {'Multi-layered sphere fit. ' ,...
'I(q) = I0*Sq*P(q; r0, sig0) + background',...
'    Sq = powI*q.^PorodExp + S(q; D, vf)',...
'    background = poly1*q.^poly2 + poly3*q + poly4',...
'CF. if you have data in absolute unit, consider use FitLee_schultzsphere2',...
'',...
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
    p = 1;
    isini = 1;
end

%% initialize fit parameter bestP
if isini
    Nf = p;
    bestP = [];
    bestP.I0 = 4.6095e-13;
    bestP.r0 = 50;
    bestP.t1 = 10;
    bestP.t2 = 10;
    bestP.t3 = 10;
    bestP.t4 = 10;
    bestP.t5 = 10;
    bestP.t6 = 10;
    bestP.t7 = 10;
    bestP.t8 = 10;
    bestP.t9 = 10;
    bestP.t10 = 10;
    bestP.ed0 = 50;
    bestP.ed1 = 10;
    bestP.ed2 = 10;
    bestP.ed3 = 10;
    bestP.ed4 = 10;
    bestP.ed5 = 10;
    bestP.ed6 = 10;
    bestP.ed7 = 10;
    bestP.ed8 = 10;
    bestP.ed9 = 10;
    bestP.ed10 = 10;
    bestP.ed_solv = 0;
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

%radius = [p.r0, p.t1, p.t2, p.t3, p.t4, p.t5, p.t6, p.t7, p.t8, p.t9, p.t10];
radius = [p.r0];
density = [p.ed0];
for i=1:10
    radius = [radius, p.(sprintf('t%i', i))];
    density = [density, p.(sprintf('ed%i', i))];
end
radius = cumsum(radius);
density = [density, p.ed_solv];

[Aq, V1] = multilayersphereamp(q, radius, density);

Pq1 = Aq(:).^2;
Sq1 = strfactor2(q, p.D, p.vf);
%Sq = p.powI*q.^p.PorodExp + Sq1;

back = p.poly1*q.^p.poly2 + p.poly3*q + p.poly4;
out = p.I0*Pq1.*Sq1 + back;
if isnan(out)
    out = ones(size(out))*1E20;
end

if nargout == 2
%     x = 0:1:radius(end);
%     rho = zeros(size(x));
%     t = x < p.coreRdivshR*p.shellR;
%     
%     rho(t) = p.core_rho;
%     t = (x >= p.coreRdivshR*p.shellR) & (x < p.shellR);
%     rho(t) = p.sh1_rho;
%     t = (x >= p.shellR) & (x < p.shellR+p.sh2thick);
%     rho(t) = p.sh2_rho;
%     t = (x >= p.shellR+p.sh2thick);
%     rho(t) = p.solvent_rho;
%     
%     figure;
%     plot(x, rho);xlabel('Radius (A)');ylabel('\rho (R)')
%     report = '';
end
