function [y, x] = multi_gauss(Rmin, Rmax, p)
% function y = multi_gauss(Rmin, Rmax, p)
% p is composed of 10 values

%R = Rmin:(Rmax-Rmin)/(prod(size(p))-1):Rmax;
%y = [];
%x = Rmin:(Rmax-Rmin)/60:Rmax;
%for i=1:prod(size(p))
%    [a, b] = gaus(x, R(i), (R(2)-R(1))/2);
%    y(:, i) = abs(p(i))*b;
%end
%y = sum(y,2);x = x';

R = Rmin:(Rmax-Rmin)/(prod(size(p))-1):Rmax;
x = Rmin:(Rmax-Rmin)/40:Rmax;
y = interp1(R, abs(p), x, 'pchip')';x = x';