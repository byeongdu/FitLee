function a = voigtarea(x, p)
% calculate the area of voigt function...
%   pa = [amp, center, sigg, sigl];
    y = abs(voigt(x, [p, 0]));
    a = trapz(x, y);
end