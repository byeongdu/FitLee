function f = saxsellipsoidCylinder(qp, qz, Rp, Rz, rho)
% cylinder with ellipsoid cross section.
% eccentricity axis is qx direction.
%R = [Rin2, Rin1, Rout];
%rho = [rho1, rho2, rho3, rho4];
% qp and qz should have a same size....

onex = ones(size(qp));
zerox = zeros(size(qz));

alpha = acos(qp./sqrt(qp.^2 + qz.^2+eps));
q = sqrt(qp.^2 + qz.^2);
%f = besseljc(q.*Re);
        
fq = zeros(size(qp));
for i=1:length(Rp)
   Re = sqrt(Rz(i)^2*sin(alpha).^2 + Rp(i)^2.*cos(alpha).^2);
   fq = fq + 2*pi*(rho(i) - rho(i+1))*Rp(i)*Rz(i)*besseljc(q.*Re);
%    fq = fq + besseljc(q.*Re);
end
f = fq;
%f = abs(fq).^2;%./(qp+eps);