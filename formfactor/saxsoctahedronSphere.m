function [F, Vwh] = saxsoctahedronSphere(qx, qy, qz, parameter)
% [F, Voch] = saxsoctahedronSphere(qx, qy, qz, [R, sphereRadius]);
% R edge length of core.
% sphereRadius :
R = parameter(1);
t = parameter(2);
rho_sh = 3.5;
rho_core = 1.4898;
[Fc, Vc] = saxsoctahedron(qx, qy, qz, R);
j = sqrt(-1);
Vs = (4*pi/3*t^3);
%xp = [0, 0, (R+t)/2, 0, -(R+t)/2, 0];
%yp = [0, 0, 0, (R+t)/2, 0, -(R+t)/2];
%zp = [sqrt(2)*(R+t), -sqrt(2)*(R+t), 0, 0, 0, 0];
p1 = [0, 0, sqrt(2)*(R+t)];
p2 = [0, 0, -sqrt(2)*(R+t)];
p3 = [(R+t)/2, (R+t)/2, 0];
p4 = [(R+t)/2, -(R+t)/2, 0];
p5 = [-(R+t)/2, (R+t)/2, 0];
p6 = [-(R+t)/2, -(R+t)/2, 0];

edges{1} = [p1;p3];
edges{2} = [p1;p4];
edges{3} = [p1;p5];
edges{4} = [p1;p6];
edges{5} = [p2;p3];
edges{6} = [p2;p4];
edges{7} = [p2;p5];
edges{8} = [p2;p6];
edges{9} = [p3;p4];
edges{10} = [p3;p5];
edges{11} = [p5;p6];
edges{12} = [p4;p6];

F = rho_core*Vc*Fc;
%qx = reshape(qx, numel(qx), 1);
%qy = reshape(qy, numel(qy), 1);
%qz = reshape(qz, numel(qz), 1);
q = sqrt(qx.^2+qy.^2+qz.^2);
p = 0.1:0.1:0.9;
for k = 1:12;
    for kk = 1:numel(p)
        [x,y,z] = equationofline(edges{k}(1,:), edges{k}(2,:), p(kk));
        F = F + rho_sh*sphereamp(q, t).*exp(-j*(qx*x + qy*y + qz*z));
    end
end
N = numel(p)*12;
Vwh = (rho_core*Vc + rho_sh*Vs*N);
F = F/Vwh;



function [x,y,z] = equationofline(p1, p2, t)
x = p1(1) + (p2(1)-p1(1))*t;
y = p1(2) + (p2(2)-p1(2))*t;
z = p1(3) + (p2(3)-p1(3))*t;