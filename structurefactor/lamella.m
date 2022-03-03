function [y, Fq, Iq] = lamella(varargin)
% y = lamella(q, thick, density, domainsize, microstrain, DW, ispowder)
q = varargin{1};
thick = varargin{2};
density = varargin{3};
domainsize = 5000;
microstrain = 0.02;
DW = 0.005;
ispowder = 0;

if numel(varargin) > 3
    domainsize = varargin{4};
end
if numel(varargin) > 3
    microstrain = varargin{5};
end
if numel(varargin) > 3
    DW = varargin{6};
end
if numel(varargin) > 3
    ispowder = varargin{7};
end
wavelength = 1; % for peakwidth calculation only.
% thick = zeros(size(relcomp));
d = sum(thick);
pos = thick;
% maxh = 5;
maxh = fix(d/(2*pi/max(q)))+2;
DW = DW*d;
Fq = zeros(1, maxh);

% for i=1:length(relcomp)
%     thick(i) = d*relcomp(i);
% end
pos(1) = 0;
j = sqrt(-1);
for i=2:length(thick)
    pos(i) = pos(i-1) + thick(i-1)/2 + thick(i)/2;
end
for h=1:maxh
    Fq(h) = 0;
    qh = 2*pi/d*h;
    for i=1:length(thick)
        Fq(h) = Fq(h) + density(i)*thick(i)*sinc(qh*thick(i)/2).*exp(-j*qh*pos(i));
    end
end
y = zeros(size(q));
y = y(:);
Iq = zeros(size(Fq));
for h = 1:maxh
    qh = 2*pi/d*h;
    Iq(h) = abs(Fq(h))^2.*exp(-qh^2*DW^2); % DW factor
    if ispowder % then, Lorentz correction.
        Iq(h) = Iq(h)./qh.^2;
    end
    w = peakwidth(qh, wavelength, domainsize, microstrain);
    t = pseudovoigt(q, [Iq(h), 2*pi/d*h, 0.00025, w]);
    y = y + t(:);
end