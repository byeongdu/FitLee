function [nr, x] = logndist99(r0, sig, numberofPoint)
%function [nr, x] = logndist99(r0, sig, numberofPoint);
% numberofPoint : # of point of x and nr...
% if sig<0.005 it will result only one.
if nargin < 3
    numberofPoint = 20 -1;
end
k=2;
if sig >= 0.005
    i=r0;
    while fix(logndistcdf2(r0, sig, i)*100)/100 < 1
        i = i+sig*r0*k;
        k = k^2
    end
    x = 0:1/1000:1;x=x*i;
%    x = 0:r0/10*(sig/0.2):r0*500*sig;
    p = logndistcdf(r0, sig, x);
    p = fix(p*100)/100;
    ninenineP = min(find(p >= 0.99));
    oneP = find(p <= 0.001);
    if isempty(oneP)
        oneP = 1;
    end
    dmin = x(oneP(1));
    dmax = x(ninenineP(1));
    x = dmin:(dmax-dmin)/numberofPoint:dmax;   %20numbers of point;
    nr = logndist(r0, sig, x);
    x=x';
else
    x=r0;nr=1;
end