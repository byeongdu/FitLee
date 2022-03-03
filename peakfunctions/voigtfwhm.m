function W = voigtfwhm(gw, lw)
% http://en.wikipedia.org/wiki/Voigt_profile
fg = 2*gw*sqrt(2*log(2));
fl = 2*lw;
W = 0.5346*fl+sqrt(0.2166*fl.^2+fg.^2);
return
% REVIEW OF SCIENTIFIC INSTRUMENTS 76, 056107 (2005)
%cs = sqrt(log(2));
%y = cs*lw./gw;
%d = (y-cs)./(y+cs);
%R = 1-0.1821*(1-d.^2)+(0.023665*exp(0.6*d)+0.00418*exp(-1.9*d)).*sin(pi*d);
%W = (y+cs).*R;