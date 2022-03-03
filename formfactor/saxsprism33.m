function [F, V] = saxsprism33(qx, qy, qz, Param)
% B. Lee book. --- something wrong with this equation....
% its center of mass is on origin.
% bottom edge of prism3 is parallel to X. 2R is the edge length. H is the
% height.
qx = qx+eps;
qy = qy+eps;
qz = qz+eps;
R = Param(1); H = Param(2);
V = sqrt(3)*R^2*H;

%t = find(abs(qy) < 0.00001);
%qy(t) = sign(qy(t))*0.00001;
%t = find(abs(qx.^2 - 3*qy.^2) < 0.00001);
%qx(t) = sqrt((0.00001-3*qy(t).^2)/3);
j = sqrt(-1);
qx = reshape(qx, size(qz));
qy = reshape(qy, size(qz));
%qx = qx+eps;
%qy = qy+eps;
F = 2*H*sqrt(3)*exp(j*qy*R/sqrt(3))./(qy.*(qx.^2 - 3*qy.^2)+eps);
F = F.*(qx.*exp(-j*qy*R*sqrt(3)) - qx.*cos(qx*R) + 3*j.*qy.*sin(qx*R));
%.*sin(q2*R/2).*exp(j*qx*R/2)-q2.*sin(q1*R/2).*exp(-j*qx*R/2));
% triangle was integrated, but H is not integrated just sinc is used.
% it should not be just sinc but H*sinc.... 
F = F.*sinc(qz*H/2);
%F = F/V;
t = isinf(F); F(t) = 1;
%t = abs(F) > 1; F(t) = 1;

t = isnan(F);
F(t) = 1;

if nargin<6
    return
end
