function [F, Vrho] = multilayersphereamp(q, radius, density)
% function [y, Vrho] = multilayersphereamp(q, radius, density)
% Vrho = V*delta_rho;

F = zeros(size(q));
Vrho = 0;

for i=1:length(radius)
    F = F + (density(i)-density(i+1))*sphereamp(q, radius(i));
	Vrho = Vrho + (density(i)-density(i+1))*radius(i).^3/3*4*pi;
end
