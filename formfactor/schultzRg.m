function [zRg, V, V2, S] = schultzRg(avgr, fwhm)
% calculate Rg from the fit parameters of schultz distribution function(r and sigma)
% assuming the spherical particles.
% [zRg, V, V2, S] = schultzRg(avgr, fwhm)
Z = (avgr./fwhm).^2 - 1;
eps = (Z+8).*(Z+7)./(Z+1).^2;
zRg =sqrt(3/5*eps).*avgr;


zp1=Z+1;
V = 4*pi/3*(zp1+1).*(zp1+2)./zp1.^2.*avgr.^3;
V2 = (4*pi/3)^2*(Z+6).*(Z+5).*(Z+4).*(Z+3).*(Z+2)./zp1.^5.*avgr.^6;
S = 4*pi*(Z+3)./(Z+2).*avgr.^2;
%sqrt((Z+8).*(Z+7))./(Z+3), %*3*sqrt(3/5)/pi
%(Z+8)^2*(Z+7)^2/((Z+6)*(Z+5)*(Z+4)*(Z+3)), %*81/50
%killwaves aaaa, bbbb
end