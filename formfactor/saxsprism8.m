function F = saxsprism8(x, y, z, R, Ra, Rb)
F = zeros(size(x)); 
qz = z;
for i=0:7
    qx = cos(pi/4*i)*x + sin(pi/4*i)*y;
    qy = -sin(pi/4*i)*x + cos(pi/4*i)*y;
    F = F + saxstriangle(qx,qy,R,pi/8);
end
if nargin > 4
    if Ra > 0
        qr = sqrt(qx.^2 + qy.^2);
%        F = F - pi*Rin^2*besseljc(qr.*Ra);
         F = F - Ra*Rb*sinc(1/pi*x*Ra).*sinc(1/pi*y*Rb);
    end
end
