function f = saxsellipsoid(qx, qy, qz, parameter)
% eccentricity axis is [1, 0, 0] or qx.
R = parameter(:,1);
e = parameter(:,2);
onex = ones(size(qx));
zerox = zeros(size(qx));
%qy = sqrt(qy.^2+qz.^2);
%alpha = acos(qx./sqrt(qx.^2 + qy.^2));
Q = [qx(:), qy(:), qz(:)];
q = sqrt(qx.^2+qy.^2+qz.^2);
alpha = angle2vect([1,0,0], Q);
Re = R*sqrt(sin(alpha).^2 + e.^2.*cos(alpha).^2);
%q = sqrt(qx.^2 + qy.^2);
f = 3*((sin(Re.*q) - Re.*q.*(cos(Re.*q))) + eps)./((Re.*q).^3 + eps);
f = f*(4*pi/3*R^3*e);