function [F, Vwh] = saxsoctahedroncoreshell(qx, qy, qz, parameter)
% [F, Voch] = saxsoctahedroncoreshell(qx, qy, qz, [R, t]);
% R edge length of core.
% t thickness of shell.
R = parameter(1);
t = parameter(2);
rho_sh = 3.5;
rho_core = 1.4898;
[Fc, Vc] = saxsoctahedron(qx, qy, qz, R);
[Fwh, Vwh] = saxsoctahedron(qx, qy, qz, R+2*t);
F = rho_sh*Vwh*Fwh - rho_sh*Vc*Fc + rho_core*Vc*Fc;
F = F/(rho_sh*Vwh-rho_sh*Vc+rho_core*Vc);