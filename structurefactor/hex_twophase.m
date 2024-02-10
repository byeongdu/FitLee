function [y, Fq, Iq] = hex_twophase(varargin)
% y = hex_twophase(q, [center_center_distance, volf_of_cylinder], delta_rho, domainsize, microstrain, DW)
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

area = a^2*sqrt(3)/2;
R = sqrt(volf*area/pi);

if numel(varargin) > 6
    R = varargin{7};
    dR = varargin{8};
    volf = pi*R^2/area;
    DW = -DW;
end

wavelength = 1; % for peakwidth calculation only.
% thick = zeros(size(relcomp));
% maxh = 5;
maxq = max(q);
maxh = fix(d/(2*pi/maxq));
maxk = maxh;
DW = DW*d;
Fq = [];
qhk = [];

%a1 = [Dsp, 0];
%a2 = [-1/2*Dsp, -sqrt(3)/2*Dsp];
%area = a^2*sqrt(3)/2;
%R = sqrt(volf*area/pi);
if volf>=0.9069
    fprintf('For hex, volf of cylinder cannot exceed 0.9069.\n');
    return
end
j = sqrt(-1);

HKL = load('bcchkl.mat');
hkl = HKL.HKL(:,1:3);
dv = HKL.HKL(:,4)/100*a;
qhk = 2*pi./dv;

%Fq = coef*cylinder_amp(qhk,0,0,R,1,[0,0,1]).*(1+...
%    exp(-j*2*pi*sum(bsxfun(@times, hkl, [1/2, 1/2, 1/2]),2)));
%Fq = Fq/(4*pi/3*R^3);

if dR ~= -1
    %Pq2 = saxscylinder_CS(qhk(:), R, [1, 0], dR);
    %Fq = Fq./sqrt(Pq2);
    Fq = ones(size(qhk(:)));
    Pq = saxscylinder_CS(q(:), R, [1, 0], dR);
else
    Pq = ones(size(q));
end

y = zeros(size(q));
y = y(:);
Iq = zeros(size(Fq));
omega = 2*pi;
factor = 1/omega*(2*pi)^2/area/10;

for h = 1:numel(qhk)
    qh = qhk(h);
    if DW>0
        Iq(h) = abs(Fq(h))^2.*exp(-qh^2*DW^2); % DW factor
    else
        Iq(h) = abs(Fq(h))^2;
    end

    %Iq(h) = multip(h)*abs(Fq(h))^2.*exp(-qh^2*DW^2); % DW factor
    if ispowder % then, Lorentz correction.
        Iq = Iq./qh;
    end
    w = peakwidth(qh, wavelength, domainsize, microstrain);
    t = pseudovoigt(q, [Iq(h), qh, 0.00025, w]);
    y = y + Iq(h)*t(:);
end
y = y*factor;
if DW<0
    y(q>2*pi/min(dv)) = 1;
    Sq = 1+(y-1).*exp(-DW^2*q.^2);
end
y = Sq.*Pq;