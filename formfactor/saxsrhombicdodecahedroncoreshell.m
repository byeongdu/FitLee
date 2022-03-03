function F = saxsrhombicdodecahedroncoreshell(qx, qy, qz, a)
% a = [core_radius, outerradius, density_core, density_shell, solvent];
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
[Fc, Vrd] = saxsrhombicdodecahedron(qx, qy, qz, a(1));
Fc = Fc*Vrd;
[Fs, Vrd] = saxsrhombicdodecahedron(qx, qy, qz, a(2));
Fs = Fs*Vrd;
F = (a(4)-a(5))*Fs+(a(3)-a(4))*Fc;