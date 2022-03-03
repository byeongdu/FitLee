function [dist, rr] = logndist9(r0, sig);
%function [nr, x] = logndist99(r0, sig, numberofPoint);
% numberofPoint : # of point of x and nr...
% if sig<0.005 it will result only one.
sig = abs(sig);
xi = r0*sig/2;
x = [];
nr = [];
sumrr = -1;
maxnr = logndist(r0, sig, r0);
while(1)
%    x = [x, xi];
    nri = logndist(r0, sig, xi);
    x = [x, xi];
    nr = [nr, nri];
    xi = xi + r0*sig/2;
    if (round(sumrr*10000000) < round(sum(nr)*10000000)) & (sum(nr)>maxnr*0.1)
        sumrr = sum(nr);
        rr = x;
        dist = nr;
    else
%        x(end) = [];
%        nr(end) = [];
        if(xi > r0)
            break
        end
    end
end
if length(rr) > 40
    rr = rr(1:fix(length(rr)/20):end);
    dist = dist(1:fix(length(dist)/20):end);
end
dist = dist/sumrr;
