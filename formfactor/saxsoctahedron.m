function [F, Voch] = saxsoctahedron(qx, qy, qz, R)
% [F, Voch] = saxsoctahedron(qx, qy, qz, R)
% Edgelength : 2*R
[F1, Vpy1] = saxspyramid(qx, qy, qz, R);
%[F2, Vpy2] = saxspyramid(qx, qy, -qz, R);
% F = F1+F2;
%Voch = Vpy1 + Vpy2;
F = 2*real(F1);
Voch = 2*Vpy1;
