function [x, funR] = gaus(x, t, sigm)
% function [x, funR] = gauss(x, t, sigm)

g = inline('1/sqrt(2*pi)/(sigm)*exp(-1/2/sigm^2*(x-t).^2)','x','t','sigm');
funR = g(x, t, sigm+eps)';
[m,n] = size(x);
if m<n
   x = x';
end