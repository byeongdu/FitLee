function [y, Fq, Iq] = hcp_twophase(varargin)
% y = bcc_twophase(q, [center_center_distance, volf_of_cylinder], delta_rho, domainsize, microstrain, DW)
q = varargin{1};
a = varargin{2}(1);
d = a;
c = sqrt(8/3)*a;
volf = varargin{2}(2);
d_rho = varargin{3};
domainsize = 5000;
microstrain = 0.02;
DW = 0.1;
ispowder = 1;

if numel(varargin) > 3
    domainsize = varargin{4};
end
if numel(varargin) > 3
    microstrain = varargin{5};
end
if numel(varargin) > 3
    DW = varargin{6};
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
vol = a^2*sqrt(3)/2*c;
R = (volf*vol/2/(4*pi/3))^(1/3);
maxpackingfraction = (pi/(3*sqrt(2)));
if volf>= maxpackingfraction
    fprintf('For HCP, volf of sphere cannot excced %0.5f.\n', maxpackingfraction);
    return
else
    fprintf('Radius of sphere is %0.3f.\n', R);
end
j = sqrt(-1);

m = 1;
coef = d_rho;

[h,k,l]=ndgrid(-maxh:maxh, -maxk:maxk, -maxl:maxl);
h = h(:);
k = k(:);
l = l(:);
dv = 1./sqrt(4/3*(h.^2+h.*k+k.^2)/a^2+l.^2/c^2);
hkl = [h,k,l];
t = h==0 & k==0 & l==0;
hkl(t,:) = [];
dv(t) = [];

[~, index, ind2] = unique(dv, 'rows', 'stable');
% multiplicity
multip = histc(ind2, 1:max(ind2));
hkl = hkl(index, :);
dv = dv(index);
qhk = 2*pi./dv;
Fq = coef*sphereamp(qhk,R).*(...
    exp(-j*2*pi*sum(bsxfun(@times, hkl, [1/3, 2/3, 1/4]),2)) + ...
    exp(-j*2*pi*sum(bsxfun(@times, hkl, [2/3, 1/3, 3/4]),2)));


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