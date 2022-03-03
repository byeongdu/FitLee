function [nr, x] = normdist99(r0, sig, numberofPoint);
%function [nr, x] = normdist99(r0, sig, numberofPoint);
% numberofPoint : # of point of x and nr...
% if sig<0.005 it will result only one.
if nargin < 3
    numberofPoint = 20 -1;
end

if sig >= 0.005
    i=r0;
    while fix(normcdf(i, r0, sig)*10000)/10000 < 1
        i = i+sig*r0;
    end
    x = 0:i/100:i;
%    x = 0:r0/10*(sig/0.2):r0*500*sig;
    p = normcdf(x, r0, sig);
    p = fix(p*1000000000)/1000000000;
    ninenineP = min(find(p >= 0.999999999));
    oneP = find(p <= 0.001);
    if isempty(oneP)
        oneP = 1;
    end
    dmin = x(oneP(1));
    dmax = x(ninenineP(1));
    x = dmin:(dmax-dmin)/numberofPoint:dmax;   %20numbers of point;
    nr = normpdf(x, r0, sig);
    x=x';
else
    x=r0;nr=1;
end