function [F, Vrd] = saxsrhombicdodecahedron(qx, qy, qz, a)
% a is the edge length;
% V = 16/9*sqrt(3)*a^3
%   = 16/9*9*R^3
%   = 16*R^3
% S = 8*sqrt(2)*a^2;
% V = 6*4/3*R^3 + (2R)^3
%   = 8*R^3+8*R^3
%   = 16*R^3
% Lazzari's paper define F(q) = Int: exp(iqr)dr instead of 
%   F(q) = Int: exp(-iqr)dr. Thus, I multiply -1 to the vector q.
% Note that in this code phase term is different from his..

%a = a/2;
j = sqrt(-1);
R = a/sqrt(3);
F1 = saxspyramid(qx, qy, qz, [R, R]);F1=F1.*exp(-j*qz*R);
%F2 = saxspyramid(qx, qy, -qz, [R, R]);F2=F2.*exp(j*qz*R);
F3 = saxspyramid(qx, qz, qy, [R, R]);F3=F3.*exp(-j*qy*R);
%F4 = saxspyramid(qx, qz, -qy, [R, R]);F4=F4.*exp(j*qy*R);
F5 = saxspyramid(qz, qy, qx, [R, R]);F5=F5.*exp(-j*qx*R);
%F6 = saxspyramid(qz, qy, -qx, [R, R]);F6=F6.*exp(j*qx*R);
%Fpy = (F1+F2+F3+F4+F5+F6);
Fpy = 2*real(F1 + F3 + F5);
Fcube = saxscube(qx,qy,qz,R);
F = Fpy + Fcube;
%F = Fpy;
Vrd = 16*R.^3;