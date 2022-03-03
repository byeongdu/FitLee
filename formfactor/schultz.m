function nr = schultz(avgr, fwhm, r)
% function nr = schultz(avgr, fwhm, r_vector)
zp1 = (avgr/fwhm)^2;
r = r(:);
if gamma(zp1) > 1E150
    [~, nr] = gaus(r, avgr, fwhm);nr = nr(:);
else
    nr = (zp1/avgr).^zp1.*r.^(zp1-1).*exp(-zp1/avgr.*r)./gamma(zp1);
end
nr(isinf(nr)) = 0;
nr(isnan(nr)) = 0;
if sum(nr) == Inf
    nr = ones(length(r), 1);
end
if isnan(sum(nr))
    nr = ones(length(r), 1);
end