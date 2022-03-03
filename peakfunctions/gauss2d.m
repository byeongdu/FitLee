function Ig=gauss2d(x, y, p)
% y=gauss2d(x, y, p)      : bivariant normal distribution function
% http://en.wikipedia.org/wiki/Multivariate_normal_distribution
% A = p(1);
% x0 = p(2);
% y0 = p(3);
% sigx = p(4);
% sigy = p(5);
% rho = p(6); % correlation between x and y.
% Area of this 2D gaussian is 1
%
% Author: B. Lee
% when p(1) =1, the area should be 1.

A = p(1);
x0 = p(2);
y0 = p(3);
sigx = p(4);
sigy = p(5);
rho = p(6); % correlation between x and y.
if sigx == 0
    Ig = ones(size(x));
    return;
end
%c = p(7);
%c1 = p(8);
%c2 = p(9);
%c3 = p(10);
%c4 = p(11);
Ig = A/(2*pi*sigx*sigy*sqrt(1-rho.^2))*exp(-1/(2*(1-rho^2))*((x-x0).^2/(sigx^2)+(y-y0).^2/(sigy^2)-2*rho*(x-x0).*(y-y0)/sigx/sigy));
%y = y+c+c1*x+c2*y+c3*x.^2+c4*y.^2;