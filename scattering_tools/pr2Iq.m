function Iq = pr2Iq(r, pr, q)
% I(q) = 4pi int pr(r) * sin(qr)/(qr) dr
% pr(r) = 1/(2pi^2) int I(q) * sin(qr)/(qr) dq

r=r(:);
pr = pr(:);
Iq = zeros(size(q));
q = q(:);
for i=1:numel(q);
    Iq(i) = trapz(r, pr .* sinc(q(i)*r));
end
Iq = Iq*(4*pi);