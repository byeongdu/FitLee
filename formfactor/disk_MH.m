function re = disk_MH(q, R)

r = 0:1:R;
size_r = length(r);
size_q = length(q);
%for i = 1:size_q
%    tt = 2*pi*r(i)*besselj(0, q(i)*r);
%    re(i) = integral(r, tt);
%end

re = (1*1*(2*besselj(1, q.*R)./(q*R))).^2;