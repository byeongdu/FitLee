function y = cylinder_amp(qx, qy, qz, R, L, cyl_Lvector)
% cylinder amplitude scattering
% y = cylinder_amp(qx, qy, qz, R, L, cyl_Lvector)
% where length of qx, qy, qz should be the same.
% cyl_Lvector is the vector parallel to the cylinder axis;
% vertical cylinder [0, 0, 1];
% parallel cylinder to x-ray [1, 0, 0];
% perpendicular cylinder to x-ray [0, 1, 0];
%cyl_Lvector = [0, 0, 1];
qvector = [qx(:), qy(:), abs(qz(:))];
cyl_Lvector = repmat(cyl_Lvector, numel(qvector)/3, 1);
alpha = angle2vect(cyl_Lvector, qvector);
q = sqrt(qx.^2+qy.^2+abs(qz).^2);q = q(:);
%centofcyl = 0;
y = 2*pi*R.^2*L*besseljc(R*q.*sin(alpha)).*sinc(L*q.*cos(alpha)/2);%.*exp(-j*qz*centofcyl);
y = reshape(y, size(qx));