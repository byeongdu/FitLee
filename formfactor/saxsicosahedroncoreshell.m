function [F, V] = saxsicosahedroncoreshell(qx, qy, qz, L)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%a = 2;
L1 = L(1);
L2 = L(2);
rho1 = L(3);
rho2 = L(4);
[F1, V1] = saxsicosahedron(qx, qy, qz, L1);
[F2, V2] = saxsicosahedron(qx, qy, qz, L2);
F = rho1*F1 - rho2*F2;
V = V1 - V2;