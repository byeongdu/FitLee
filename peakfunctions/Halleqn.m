function Intensity = Halleqn(vari, s)
% Intensity = Halleqn(vari, s)
% vari = [N, deltaRo, NuY, NuZ, Y, Z]

N = vari(1);
deltaRo = vari(2);
NuY = vari(3);
NuZ = vari(4);
Y = vari(5);
Z = vari(6);

A = 1./(1+(2*pi*s*NuY*Y).^2);
B = 1./(1+(2*pi*s*NuZ*Z).^2);

epsiY = Y*(1-2*NuY);
epsiZ = Z*(1-2*NuZ);

phiY = (2*pi*s*epsiY + acos((1-(2*pi*s*NuY*Y).^2)./(1+(2*pi*s*NuY*Y).^2)));
phiZ = (2*pi*s*epsiZ + acos((1-(2*pi*s*NuZ*Z).^2)./(1+(2*pi*s*NuZ*Z).^2)));

X = phiY + phiZ;

F = 1+(A.^2).*(B.^2) - 2*A.*B.*cos(X);
C = N*(1-(A.^2).*(B.^2) - A.*(1-B.^2).*cos(phiZ) - B.*(1-A.^2).*cos(phiY));
G = 1-(A.^N).*(B.^N).*cos(N.*X);
D = B.*((1-A.^2).*(1-(A.^2).*(B.^2)).*sin(X).*sin(phiZ) + ((1+(A.^2).*(B.^2)).*cos(X) - 2*A.*B) +((1+A.^2).*cos(phiZ)-2*A))./G;
E = (A.^N).*(B.^(N+1)).*sin(N.*X).*((1-(A.^2).*(B.^2)).*((1+(A.^2).*cos(phiZ) - 2*A)).*sin(X) -(1-(A.^2)).*((1+A.^2.*(B.^2)).*cos(X) - 2*A.*B).*sin(phiZ));

Intensity = (C./F + (D+E)./F.^2)./(pi*s).^2/deltaRo^2;