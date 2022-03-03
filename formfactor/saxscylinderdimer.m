function [y, Vcyl2] = saxscylinderdimer(qx, qy, qz, p)
% function [y, Vcyl2] = saxscylinderdimer(qx, qy, qz, p)
% p = [R, L, tip-to-tip gap, dihedralangle]
RL = p(1:2);
F1 = saxscylinder2(qx, qy, qz, RL);
theta = p(4);
R = rotate_around_vector([1, 0, 0], theta);
Q = [qy(:), qz(:)]*R(2:end, 2:end);
qyn = Q(:,1); qzn = Q(:,2);
qyn = reshape(qyn, size(qx));
qzn = reshape(qzn, size(qx));
[F2, Vcyl2] = saxscylinder2(qx, qyn, qzn, RL);
Vcyl2 = Vcyl2*2;
%yshift = p(2)/2*sin(theta*pi/180);
%zshift = p(2)/2*(1-cos(theta*pi/180));
xshift = p(3)+2*p(1);
%F2 = F2.*exp(-sqrt(-1)*qy*yshift).*exp(-sqrt(-1)*qz*zshift).*exp(-sqrt(-1)*qx*xshift);
F2 = F2.*exp(-sqrt(-1)*qx*xshift);
y = F1+F2;

