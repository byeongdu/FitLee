function [nr, x] = schultzdist99(r0, sig, numberofPoint)
%function [nr, x] = schultzdist99(r0, sig, numberofPoint);
% numberofPoint : # of point of x and nr...
% if sig<0.001 it will result only one.
if nargin < 3
    numberofPoint = 20 -1;
end
if (numel(r0) == 1)
    sig = sig(1);
end
if r0 ==0
    nr = 1;
    x = 0;
    return
end
if sig/r0 < 0.001
    nr = 1;
    x = r0;
    return
end
if sig >= r0/100
    i=r0;
    while fix(schultzdistcdf(r0, sig, i)*100)/100 < 0.99
        i = i+sig/r0*10;
        if i>1E4
            i=1E4;
            break
        end
    end
    x = 0:i/100:i;
    p = schultzdistcdf(r0, sig, x);
    p = fix(p*100)/100;
    t = find(p >= 0.99);ninenineP = min(t);
    t = find(p <= 0.01);oneP = max(t);
    if isempty(oneP)
        oneP = 1;
    end
    dmin = x(oneP(1));
    dmax = x(end);
    if ninenineP > 0
        dmax = x(ninenineP(1));
    end
    x = dmin:(dmax-dmin)/(numberofPoint-1):dmax;   %20numbers of point;
    nr = schultzdist(x, r0, sig);
else
    x=r0;nr=1;
end