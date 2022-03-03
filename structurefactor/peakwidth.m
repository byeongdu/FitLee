function w = peakwidth(q, wavelength, dsize, mstrain)
% w = peakwidth(q, wavelength, dsize, mstrain)
% Calculated the width of a diffraction that is Lorenzian,
% when domain size and micro-strain are given. 
% See Senesi's paper and C. Weidenthaler, nanoscale, 2011, 3, 792;


%ang = q2angle(q, wavelength)/2;
kv = 0.9;
%w = (k/dsize*wavelength + 4*mstrain*sin(ang*pi/180))/cos(ang*pi/180);
% beta = 4*(pi*k/L+epis*q)*lambda/sqrt(16*pi^2-(q*lambda).^2)
% Integral breadth
w = 4*(pi*kv/dsize+mstrain*q)*wavelength/sqrt(16*pi^2-(q*wavelength).^2);
% FWHM of Lorentzian = (pi/2)*FWHM
w = w*pi/2;
end
