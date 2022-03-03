function [ysum, V2sum] = formfactorscale(q, Iq, p, newq, newparameter, dr, vlmscale)
% ysum = formfactorscale(q, Iq, p, newq, newparameter, dr, vlmscale)
% p : original size to calculate Iq vs q
% newq: new q for ysum
% newparameter: r0 (or peak of new size distribution)
% dr: width of the size distribution function
% scale the form factor calculated for a particle with size of p to that
% for a particle with the size of newparamter.
% This will be useful for calculating SAXS fo randomly oriented cubes.
% When the form factor of I(0) =1, you need to multiply volume^2 to form
% factor.

if nargin < 7
    vlmscale = 0;
end

oq = q*p;
%nq = newq*newparameter;
%newIq = interp1(oq, Iq, nq);
dx = 2*pi/max(newq);
NumX = (fix(dr*9/dx)+1)*3;
if NumX < 30
    NumX = 30;
end
if NumX > 150
    NumX = 150;
    fprintf('Number points for the distribution calculation is reduced to %i\n', NumX);
end
if (dr>0)
    x = linspace(max([0, newparameter-dr*8]), newparameter*3, NumX);
    x(x < 1) = [];

    %[x, nr] = gaus(x, newparameter, dr);
    nr = schultz(newparameter, dr, x);
    dx = (x(2)-x(1));
else
    x = newparameter;
    dx = 1;
    nr = 1;
end
V2sum = 0;
ysum = 0;

Iq_V2normalized = Iq/p^6;
%yy = [];
for m=1:numel(nr)
    %qn = newq*x(m);
    %qn = oq/x(m);
    oq2 = oq/x(m);
    yn = interp1(oq2, Iq_V2normalized, newq);
%    yn = spline(oq2, Iq_V2normalized, newq);
    %yy = [yy, yn];
    %yn = interp1(oq, Iq, qn);
    if sum(isnan(yn))>1
        fprintf('Warning! Not enough high q data for formfactorscale.m\n');
        yn(isnan(yn)) = 0;
    end
    if vlmscale
        yn = yn*(x(m)/p)^6; % volume ratio. assuming 3D particles.
    end
    ysum = ysum + nr(m)*yn*dx;
    V2sum = V2sum + nr(m)*(x(m)/p)^6*dx;
end
ysum = ysum*p^6;%*V2sum;

%ysum = ysum/sum(nr.^2);
