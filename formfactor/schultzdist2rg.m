function Rg = schultzdist2rg(varargin)
% [nV, nV2] = dist2Rg(x, nDist)
% from number distribution, it calculate Rg, volume distribution, and
% volume avergae particle radius, V^2 average particle radius

r1 = varargin{1};
fwhm1 = varargin{2};
I01 = varargin{3};
r2 = varargin{4};
fwhm2 = varargin{5};
I02 = varargin{6};
x = 0.1:0.1:500;
x = x(:);

Rg = zeros(size(r1));
for i=1:numel(r1)
    dist = schultzdist(x, r1(i), fwhm1(i))*I01(i) + schultzdist(x, r2(i), fwhm2(i))*I02(i);
    dist = dist(:)';
    Rg(i) = sqrt(3/5*(dist*(x.^8))/(dist*(x.^6)));  % sphere
end