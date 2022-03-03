function pr = Iq2pr(q, Iq, r)
% I(q) = 4pi int pr(r) * sin(qr)/(qr) dr
% pr(r) = r^2/(2pi^2) int I(q) q^2 * sin(qr)/(qr) dq

r=r(:);
%pr = pr(:);
gamma = zeros(size(r));
Iq = Iq(:);
q = q(:);
for i=1:numel(r);
    gamma(i) = trapz(q, q.^2.*Iq.*sinc(q*r(i)));
end
pr = gamma/(2*pi^2).*r.^2;