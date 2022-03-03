function [Intensity, int1, int2] = lamsaxs(lc, la, N, g, q)

Fn = sin(lc*q/2)./(q/2);

L = la + lc;

Fntotal = (sin((lc*N+la*(N-1))*q/2)./(q/2)).^2;

rho = exp(-1.*g^2/2*L^2.*q.^2);
Zn = 0;

for m = 1:N
   Z = (1-m/N).*rho.^m.*cos(m.*q.*L);
   Zn = Zn + Z;
end

Zn = 1+2*Zn;

R  = normrnd(lc, g, N,1);
Fnall = Fn;

for i=1:N
   Fnall(i, :) = sin(R(i).*q/2)./(q/2);
end

F1 = mean(abs(Fnall).^2,1);
F2 = mean(Fnall, 1).^2;

int1 = Fn.^2.*Zn;
int2 = N.*F1 - N.*F2;
Intensity = Fn.^2.*Zn + N.*F1 - N.*F2;