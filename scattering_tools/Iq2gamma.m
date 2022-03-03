function gm = Iq2gamma(q, Iq, r)
% I(q) = 4pi int gamma(r) r^2 * sin(qr)/(qr) dr
% gamma(r) = 1/(4pi) int I(q)q^2 * sin(qr)/(qr) dq

dq = q(2)-q(1);
q=q(:);
Iq = Iq(:);
gm = zeros(size(r));
r = r(:);
for i=1:numel(r);
    gm(i) = sum(Iq.*q.^2 .* sinc(q*r(i)));
end
gm = gm/(4*pi)*dq;