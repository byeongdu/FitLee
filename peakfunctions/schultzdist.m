function y = schultzdist(x, avgr, sig)
% avgr : average 
% sig : standard deviation in a gaussian distribution(== root mean square
% deviation from a average)
%
% usage : y = schultzdist(x, avgr, sig);
z = (avgr/sig).^2-1;
y = ((z+1)/avgr).^(z+1)*x.^z.*exp(-(z+1)/avgr*x)./gamma(z+1);
%zz = (avgr/sig).^2;
%y = (((zz)^(zz))/avgr).*exp(-gammaln(zz)).*exp(-(zz)*x/avgr).*(x/avgr).^(zz-1);
t = isnan(y) + isinf(y);
t = find(t>0);
if ~isempty(t)
    y = 1/(sqrt(pi*2)*sig)*exp(-1/2*((x-avgr)/sig).^2);
end