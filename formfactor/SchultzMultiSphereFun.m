function y = SchultzMultiSphereFun(varargin)
% y=SchultzMultiSphereFun(varargin)
% analytical equation
%   without any varargin, the function prints default parameters.
% y=SchultzMultiSphereFun(q, r, fwhm_percent, n), where r could be an array.
% if r is an array, there should be n.
% n is the number of particles whose radius is r.
% P(q) = x_i / sum(n) * |F(q, r_i)|^2
q = varargin{1};
r = varargin{2};
fwhm = varargin{3};
if numel(varargin) > 3
    n = varargin{4};
else
    n = 1;
end
y = 0;
sumn = sum(n);
for i=1:numel(r)
    [P, v] = SchultzsphereFun(q, r(i), r(i)*fwhm);
    y = y + n(i)/sumn*P*v^2;
end
    