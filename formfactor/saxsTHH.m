function [F, V] = saxsTHH(qx, qy, qz, parameter)
% [F, V] = saxsTHH(qx, qy, qz, [R, L])
% scattering from THH
% The center of a mass of the THH is at (0,0,0)
% The half long axis (L) of the hexahedron section is along z 
% The half edge of two short axes (R and R) of the hexahedron section of THH
% are pointing x and y.
qx = round(qx*100000)/100000;
qy = round(qy*100000)/100000;
qz = round(qz*100000)/100000;
if numel(parameter) == 2
    R = parameter(1);
    L = parameter(2);
elseif numel(parameter) == 1
    R = parameter(1);
    L = R;
end

tanalpha = 3/7;
H = tanalpha * R;
i = sqrt(-1);
% rotation matrix
% M = [1,0,0;0,0,1;0,-1,0] and the wedge of anisopyramid points -y
% M = [1,0,0;0,0,-1;0,1,0] and the wedge of anisopyramid points y
% M = [0,0,-1;0,1,0;1,0,0]*[0,1,0;-1,0,0;0,0,-1] and the wedge of anisopyramid points x
% M = [0,0,1;0,1,0;-1,0,0]*[0,1,0;-1,0,0;0,0,-1] and the wedge of anisopyramid points -x
% Qnew = inv(M)*Qold
if L~=R;
    [F1, Vanpy] = saxsanisopyramid(qx, -qz, qy, [R, L, H]);F1 = F1.*exp(-i*qy.*R);
    F2 = saxsanisopyramid(qx, qz, -qy, [R, L, H]); F2 = F2.*exp(i*qy.*R);
    F3 = saxsanisopyramid(-qy, qz, -qx, [R, L, H]); F3 = F3.*exp(i*qx.*R);
    F4 = saxsanisopyramid(-qy, -qz, qx, [R, L, H]); F4 = F4.*exp(-i*qx.*R);
else
    [F1, Vanpy] = saxspyramid(qx, -qz, qy, [R, H]);F1 = F1.*exp(-i*qy.*R);
    F2 = saxspyramid(qx, qz, -qy, [R, H]);F2 = F2.*exp(i*qy.*R);
    F3 = saxspyramid(-qy, qz, -qx, [R, H]);F3 = F3.*exp(i*qx.*R);
    F4 = saxspyramid(-qy, -qz, qx, [R, H]);F4 = F4.*exp(-i*qx.*R);
end
F = F1+F2+F3+F4;
%F = F1;
[F1, Vpy] = saxspyramid(qx, qy, qz, [R, H]);
F1 = F1.*exp(-i*qz*L);
F2 = saxspyramid(qx, qy, -qz, [R, H]); F2 = F2.*exp(i*qz*L);
F = F + F1+F2;

[Fc, Vcube] = saxscube(qx, qy, qz, [R, R, L]);
F = F + Fc;
V = Vanpy*4+Vpy*2+Vcube;