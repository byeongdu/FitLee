function [F, V] = saxsassymcoreshell(qx, qy, qz, parameter)
% saxsassymcoreshell(qx, qy, qz, [R1, R2, Pz, contrast1, contrast2)

if numel(parameter) == 5
    R1 = parameter(1);
    R2 = parameter(2);
    Pz = parameter(3); % position of the R2 sphere along the z direction
    con1 = parameter(4); % (density - solvent)
    con2 = parameter(5); % (density - solvent)
elseif numel(parameter) == 3
    R1 = parameter(1);
    R2 = parameter(2);
    Pz = parameter(3); % position of the R2 sphere along the z direction
    con1 = 1;
    con2 = 1;
else
    error('parameter should contain 5 or 3 values, R1, R2, Pz, contrast1, contrast2')
end
if R2+Pz > R1
    error('R2 + Pz should not be larger than R1')
end
j = sqrt(-1);
q = sqrt(qx.^2+qy.^2+qz.^2);
F = sphereamp(q, R1);
F2 = sphereamp(q, R2).*exp(-j*qz*Pz);

F = con1*F - con2*F2;
V = 4*pi/3*(R1^3-R2^3);