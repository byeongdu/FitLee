function [y, Fq, Iq] = bcc_twophase(varargin)
% y = bcc_twophase(q, [center_center_distance, volf_of_cylinder], delta_rho, domainsize, microstrain, DW)
q = varargin{1};
a = varargin{2}(1);
d = a;
volf = -1;
if numel(varargin{2})==2
    volf = varargin{2}(2);
end

d_rho = varargin{3};
domainsize = 5000;
microstrain = 0.02;
DW = 0.1;
ispowder = 1;
dR = -1;

if numel(varargin) > 3
    domainsize = varargin{4};
end
if numel(varargin) > 3
    microstrain = varargin{5};
end
if numel(varargin) > 3
    DW = varargin{6};
end

vol = a^3;
R = (volf*vol/2/(4*pi/3))^(1/3);

if numel(varargin) > 6
    R = varargin{7};
    dR = varargin{8};
    volf = 3*pi/3*R^3/a^3;
    DW = -DW;
end

wavelength = 1; % for peakwidth calculation only.
% thick = zeros(size(relcomp));
% maxh = 5;
maxq = max(q);
maxh = fix(d/(2*pi/maxq));
maxk = maxh;
maxl = maxh;
DW = DW*d;
Fq = [];
qhk = [];

%a1 = [Dsp, 0];
%a2 = [-1/2*Dsp, -sqrt(3)/2*Dsp];
if volf>=0.6802
    fprintf('For BCC, volf of sphere cannot excced 0.6802.\n');
    return
else
    fprintf('Radius of sphere is %0.3f.\n', R);
end
j = sqrt(-1);

m = 1;
coef = d_rho;
HKL = load('bcchkl.mat');
hkl = HKL.HKL(:,1:3);
dv = HKL.HKL(:,4)/100*a;
% [h,k,l]=ndgrid(-maxh:maxh, -maxk:maxk, -maxl:maxl);
% h = h(:);
% k = k(:);
% l = l(:);
% dv = a./sqrt(h.^2+k.^2+l.^2);
% hkl = [h,k,l];
% t = h==0 & k==0 & l==0;
% hkl(t,:) = [];
% dv(t) = [];

%[~, index, ind2] = unique(dv, 'rows', 'stable');
%% multiplicity
%multip = histc(ind2, 1:max(ind2));
%hkl = hkl(index, :);
%dv = dv(index);
qhk = 2*pi./dv;
Fq = coef*sphereamp(qhk,R).*(1+...
    exp(-j*2*pi*sum(bsxfun(@times, hkl, [1/2, 1/2, 1/2]),2)));
Fq = Fq/(4*pi/3*R^3);
if dR ~= -1
    Pq2 = SchultzsphereFun(qhk, R, dR*R);
    Fq = Fq./sqrt(Pq2);
    Pq = SchultzsphereFun(q, R, dR*R);
else
    Pq = ones(size(q));
end

y = zeros(size(q));
y = y(:);
Iq = zeros(size(Fq));
omega = 4*pi;
cVol = a^3;
factor = 1/omega*(2*pi)^3/cVol*1000;
for h = 1:numel(qhk)
    qh = qhk(h);
    if DW>0
        Iq(h) = abs(Fq(h))^2.*exp(-qh^2*DW^2); % DW factor
    else
        Iq(h) = abs(Fq(h))^2;
    end

    %Iq(h) = multip(h)*abs(Fq(h))^2.*exp(-qh^2*DW^2); % DW factor
    if ispowder % then, Lorentz correction.
        Iq = Iq./qh.^2;
    end
    w = peakwidth(qh, wavelength, domainsize, microstrain);
    t = pseudovoigt(q, [1, qh, 0.00025, w]);
    y = y + Iq(h)*t(:);
end
y = y*factor;
if DW<0
    y(q>2*pi/min(dv)) = 1;
    Sq = 1+(y-1).*exp(-DW^2*q.^2);
end
y = Sq.*Pq;