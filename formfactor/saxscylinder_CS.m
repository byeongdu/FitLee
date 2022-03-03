function Iq = saxscylinder_CS(qp, R, density, sigma)
% This is SAXS for radial direction of a Core-Shell cylinder
% R = [Rcore, Rshell];
% density = [edensity_core, edensity_shell, edensity_solvent];
% see also saxscylinder.m
% see also cylinder_amp.m, cylinder_vert_type.m
% p = [R, L]

%H = p(:,2);
%Vcyl = pi*R.^2;
if nargin < 4
    sigma = 0;
end
if sigma == 0
    F = zeros(size(qp));
    for i=1:numel(R)
        Ft = saxscylinder3(qp, R(i));
        F = F + (density(i+1)-density(i))*Ft;
    end
    Iq = abs(F).^2;
    return
end

numpnt = 41;
maxR = max(R);
x = linspace(sqrt(maxR/10), sqrt(maxR + maxR*sigma*9), numpnt);x = x(:).^2;
nr = schultzdist(x, maxR, maxR*sigma);
nr = nr(:);
Iq = zeros(size(qp));
for i=1:numpnt
    y = saxscylinder_CS(qp, R/maxR*x(i), density);
    Iq = Iq + nr(i)*y*(x(2)-x(1));
end