function [F, Vpy] = saxspyramid0(qx, qy, qz, RH)
% function [F, Vpy] = saxspyramid(qx, qy, qz, RH)
L = RH(1)*2;
R = L*sqrt(2)/2;
H = RH(2);
F = zeros(size(qx));
Vpy = L.^2.*H/3;
for i=0:3
    M = [cos(2*pi/4*i),sin(2*pi/4*i),0;-sin(2*pi/4*i),cos(2*pi/4*i),0;0,0,1];
    q1 = M(1,1)*qx + M(1,2)*qy;
    q2 = M(2,1)*qx + M(2,2)*qy;
    F1 = saxsoctant(q1, q2, qz, [R, L, H]);

    F = F+F1;
end
%F = F;