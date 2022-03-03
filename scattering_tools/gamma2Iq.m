function Iq = gamma2Iq(r, gamma, q)
% I(q) = 4pi int gamma(r) r^2 * sin(qr)/(qr) dr
% gamma(r) = 1/(2pi^2) int I(q)q^2 * sin(qr)/(qr) dq

dr = r(2)-r(1);
r=r(:);
gamma = gamma(:);
Iq = zeros(size(q));
q = q(:);
for i=1:numel(q);
    Iq(i) = sum(gamma.*r.^2 .* sinc(q(i)*r));
end
Iq = Iq*(4*pi)*dr;