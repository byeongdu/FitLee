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
    volf = pi*R^2/area^3;
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
else
    fprintf('Radius of cylinder is %0.3f.\n', R);
end
j = sqrt(-1);

m = 1;
coef = d_rho*pi*R^2;
for h=-maxh:maxh
    for k=-maxk:maxk
        dv = d./sqrt(4/3*(h.^2+h.*k+k.^2));
        qv = 2*pi/dv;
        if h==0 & k==0
            continue;
        end
        if qv > maxq
            continue;
        end
        qhk(m) = qv;
        Fq(m) = 0;
%        for i=1:length(thick)
            Fq(m) = Fq(m) + coef*besseljc(qv*R);
%        end
        m = m + 1;
    end
end

[qhk, ind, ic1] = unique(qhk);
multip = zeros(size(qhk));
for i=1:numel(qhk)
    multip(i) = sum(ic1==i);
end

if dR ~= -1
    Pq2 = saxscylinder_CS(qhk, R, [1, 0], dR);
    Fq = Fq./sqrt(Pq2);
    Pq = saxscylinder_CS(q(:), R, [1, 0], dR);
else
    Pq = ones(size(q));
end

Fq = Fq(ind);
y = zeros(size(q));
y = y(:);
Iq = zeros(size(Fq));
for h = 1:numel(qhk)
    qh = qhk(h);
    Iq(h) = multip(h)*abs(Fq(h))^2.*exp(-qh^2*DW^2); % DW factor
    if ispowder % then, Lorentz correction.
        Iq = Iq./qh.^2;
    end
    w = peakwidth(qh, wavelength, domainsize, microstrain);
    t = pseudovoigt(q, [Iq(h), qh, 0.00025, w]);
    y = y + t(:);
end
y = y*Pq;