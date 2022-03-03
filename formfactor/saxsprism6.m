function F = saxsprism6(x, y, z, parameters)
%qx = x*sqrt(3)/2+1/2*y;
%qy = -1/2*x + sqrt(3)/2*y;
%qz = z;
%F = saxsprism3(qx, qy, qz, R, H) ...
%    + saxsprism3(-qx, -qy, qz, R, H) ...
%    + saxsprism3(1/2*qx-sqrt(3)/2*qy, sqrt(3)/2*qx+1/2*qy, qz, R, H) ...
%    + saxsprism3(-1/2*qx+sqrt(3)/2*qy, -sqrt(3)/2*qx-1/2*qy, qz, R, H) ...
%    + saxsprism3(-1/2*qx-sqrt(3)/2*qy, sqrt(3)/2*qx-1/2*qy, qz, R, H) ...
%    + saxsprism3(1/2*qx+sqrt(3)/2*qy, -sqrt(3)/2*qx+1/2*qy, qz, R, H);
R = parameters(1);
H = parameters(2);
F = zeros(size(x));
qz = z;
for i=0:5
    qx = cos(pi/3*i)*x + sin(pi/3*i)*y;
    qy = -sin(pi/3*i)*x + cos(pi/3*i)*y;
%    F = F + saxsprism3(qx, qy, qz, [R, H]).*exp(-j*x*R).*exp(-j*y*R/sqrt(3));
%   shift saxsprism3 down along -y by 2*R/sqrt(3).. changed 2010. 5. 6.
    F = F + saxsprism3(qx, qy, qz, [R, H]).*exp(j*y*2*R/sqrt(3));
end