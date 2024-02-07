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

[h,k,l]=ndgrid(-maxh:maxh, -maxk:maxk, -maxl:maxl);
h = h(:);
k = k(:);
l = l(:);
dv = a./sqrt(h.^2+k.^2+l.^2);
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
Fq = coef*sphereamp(qhk,R).*(1+...
    exp(-j*2*pi*sum(bsxfun(@times, hkl, [1/2, 1/2, 1/2]),2)));

if dR ~= -1
    Pq2 = SchultzsphereFun(qhk, R, dR*R);
    Fq = Fq./sqrt(Pq2);
    Pq = SchultzsphereFun(q, R, dR*R);
else
    Pq = ones(size(q));
end

% for h=-maxh:maxh
%     for k=-maxk:maxk
%         for l=-maxl:maxl
%             dv = d./sqrt(h^2+k^2+l^2);
%             qv = 2*pi/dv;
%             if h==0 & k==0 & l==0
%                 continue;
%             end
%             if qv > maxq
%                 continue;
%             end
%             qhk(m) = qv;
%             Fq(m) = coef*sphereamp(qv,R)*(1 + exp(-j*2*pi*dot([h, k, l], [1/2, 1/2, 1/2])));
%             m = m + 1;
%         end
%     end
% end
% 
% [qhk, ind, ic1] = unique(qhk);
% multip = zeros(size(qhk));
% for i=1:numel(qhk)
%     multip(i) = sum(ic1==i);
% end
% Fq = Fq(ind);
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
y = y.*Pq;