function [F, V] = saxsprism3(qx, qy, qz, Param)
% J. appl. cyrst (2002), 35, 406, Lazzari
% its center of mass is on origin.
% bottom edge of prism3 is parallel to X. 2R is the edge length. H is the
% height.
% Lazzari's paper define F(q) = Int: exp(iqr)dr instead of 
%   F(q) = Int: exp(-iqr)dr. Thus, I multiply -1 to the vector q.
%   like below...
%   qx = -qx(:);
%   qy = -qy(:);
%   qz = -qz;
% But, I found that there was typo in the eqatuion.
%   So, I rederived the equation with F(q) = Int: exp(-iqr)dr definition.
% 11/25/2012
% B. Lee
R = Param(1); H = Param(2);
V = sqrt(3)*R^2*H;

if abs(qx) < 0.00001
    qx = 0.001*ones(size(qx));
end
t = find(abs(qx) < 0.00001);
qx(t) = sign(qx(t))*0.00001;
t = find(abs(qx.^2 - 3*qy.^2) < 0.00001);
qy(t) = sqrt((qx(t).^2 - 0.00001)/3);
j = sqrt(-1);
qx = reshape(qx, size(qz));
qy = reshape(qy, size(qz));
%qx = qx+eps;
%qy = qy+eps;
%F = 2*H*sqrt(3)*exp(-j*qy*R/sqrt(3))./(qx.*(qx.^2 - 3*qy.^2));
%F = F.*(qx.*exp(j*qy*R*sqrt(3)) - qx.*cos(qx*R) - j*sqrt(3).*qy.*sin(qx*R));

F = -2*H*sqrt(3)*exp(j*qy*R/sqrt(3))./(qx.*(qx.^2 - 3*qy.^2));
F = F.*(qx.*exp(-j*qy*R*sqrt(3)) - qx.*cos(qx*R) + j*sqrt(3).*qy.*sin(qx*R));
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