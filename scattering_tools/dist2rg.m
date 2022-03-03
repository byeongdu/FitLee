function [RgSphere, nV, nV2, meanV] = dist2rg(x, nDist)
% [nV, nV2] = dist2Rg(x, nDist)
% from number distribution, it calculate Rg, volume distribution, and
% volume avergae particle radius, V^2 average particle radius

if length(x) == 1
    [nDist, x] = logndist9(x, nDist);
end
[nx, ny] = size(nDist);
if ny < nx
    nDist = nDist';     % make nDist as column vector
end
nDist = nDist/sum(nDist);    % make sum of nDist as 1
[nx, ny] = size(x);
if ny > nx
    x = x';         % make x as row vector
end
r_num_avg = nDist*x;
meanV = sum(nDist.*(x'.^3))./sum(nDist);
nV = (nDist.*(x'.^3))./sum(nDist.*(x'.^3));
r_vol_avg = (nDist*(x.^3)).^(1/3);
nV2 = (nDist.*(x'.^6))./sum(nDist.*(x'.^6));
r_vol2_avg = (nDist*(x.^6)).^(1/6);
RgPlate = sqrt(1/2*(nDist*(x.^8))/(nDist*(x.^6))); % cylinder
RgPlate2 = sqrt(1/2*(nDist*(x.^6))/(nDist*(x.^4))); % cylinder
RgSphere = sqrt(3/5*(nDist*(x.^8))/(nDist*(x.^6)));  % sphere
sprintf('Rg for plate and sphere is %0.5g %0.5g %0.5g', RgPlate, RgPlate2, RgSphere)
%sprintf('Rg for plate and sphere is %0.5g', RgSphere)
sprintf('number averaged Radius is %0.5g', r_num_avg);
sprintf('volume averaged Radius is %0.5g', r_vol_avg);
sprintf('volume square averaged Radius is %0.5g', r_vol2_avg);