function Fq = rho2Fq(r, rho, q)
% For a spherically symmetric object either by geometry or random
% orientation,
% F(q) = int_0^inf 4*pi*r^2*rho(r)*sinc(qr) dr
%, where rho(r) is the radial scattering length density.
% 
% Byeongdu Lee, 3/2/2016
% Made for BOXmodel.m and boxPq.m

% check whether data is equally spaced or not.
z=linspace(r(1),r(end),length(r));z=z';
r=r(:);
if ~all(abs(r-z)<=abs(z.*eps))
    r2 = z;
    rho2 = interp1(r, rho, r2);
    r = r2(:);
    rho = rho2(:);
end

for i=1:numel(q)
    Fq(i) = trapz(r, 4*pi*r.^2.*rho.*sinc(q(i)*r));
end