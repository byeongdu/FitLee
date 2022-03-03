function F = saxsprism8o(x, y, z, R, Rin)
F = zeros(size(x));
dR = R*0.1;
if (dR ~= 0)
     [gDistR, Rout] = schultzdist99(R, dR, 10);gDistR=gDistR*(Rout(2)-Rout(1));
end
if nargin < 5
    Rin = 0;
end
for i=1:length(Rout)
%    qx = cos(pi/36*i)*x + sin(pi/36*i)*y;
%    qy = -sin(pi/36*i)*x + cos(pi/36*i)*y;
    R = Rout(i);
    F = F + gDistR(i)*abs(saxsprism8(x, y, z, R, Rin)).^2;
end
