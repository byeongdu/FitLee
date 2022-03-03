function F = saxsprism32(qx, qy, qz, R, H)
% J. appl. cyrst (2002), 35, 406
% Lazzari
qx = qx+eps;
qy = qy + eps;
%F = (qx.*exp(i.*(-qx.*R+qy*sqrt(3).*R-qy*sqrt(3)))-exp(-i.*R.*qx).*qx+exp(-i.*R.*qx).*qy*sqrt(3)-exp(-i.*qy*sqrt(3)).*qx+qx-qy*sqrt(3))./(-qx+qy*sqrt(3))./qy./qx-(exp(-i.*qy*sqrt(3)).*qx-qx-qy*sqrt(3)-qx.*exp(i.*(qx.*R-qy*sqrt(3)+qy*sqrt(3).*R))+exp(i.*R.*qx).*qx+exp(i.*R.*qx).*qy*sqrt(3))./(qx+qy*sqrt(3))./qx./qy;
F = (-qx.*exp(-i*R*(qx+qy*sqrt(3)))+exp(-i*R*qx).*qx+exp(-i*R*qx).*qy*sqrt(3)-qy*sqrt(3))./(qx+qy*sqrt(3))./qy./qx+(qx.*exp(2*i*(-qx*R+qy*sqrt(3)*R-qy*sqrt(3)))-exp(-2*i*R.*qx).*qx+exp(-2*i*R.*qx).*qy*sqrt(3)-exp(-2*i*qy*sqrt(3)).*qx+qx-qy*sqrt(3))./(-qx+qy*sqrt(3))./qy./qx;
 
%F = F.*exp(-i*R*qx);