function [nr, x] = logndist(r0, sig, x)
% [nr, x] = logndist(r0, sig, x);
% lognormal function
x(find(x==0)) = eps;
nr = 1/(sqrt(2*pi).*r0.*sig.*exp(0.5*sig.^2)).*exp(-1*log(x./r0).^2./(2*sig.^2));
% cdf of this function looks not analytically calculable...
%mu = log(r0);
%nr = 1./(sig*sqrt(2*pi)*x).*exp(-1*(log(x) - mu).^2./(2*sig.^2)); % matlab
%implemented lognormal distribution...
%x=x';nr=nr';