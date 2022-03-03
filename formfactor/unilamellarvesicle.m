function Fq2 = unilamellarvesicle(q, R0, eps, rho, sig)
% function Fq2 = unilamellarvesicle(q, R0, eps, rho, sig)
lenshell = length(eps);
y = zeros(size(q));
for k=1:lenshell
    for j=1:lenshell
        y = y + (R0+eps(j))*(R0+eps(k))*rho(j)*rho(k)*sig(j)*sig(k)*exp(-q.^2*(sig(k)^2+sig(j)^2/2)).*cos((q*(eps(j)-eps(k))));
    end
end
Fq2 = y./q.^2;