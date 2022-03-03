function test = OZ(vari, q)
% OZ(vari, q)
% I0 = vari(1);
% xi = vari(2);
% C = vari(3) : compressibility effect.

I0 = vari(1);
xi = vari(2);
C = vari(3);
test = vari(1)./(1+xi.^2.*q.^2)+C;