function Fq = parallelpiped(qx, qy, qz, W, L, H)
% Fq = parallelpiped(qx, qy, qz, W, L, H, beta)
% beta is the angle between axis of cylinder and qy
qx = real(qx);
qy = real(qy);
qz = real(qz);
V = W*L*H;
%size(qx), size(qy), size(qz)
Fq = V*sinc(qx*L/2).*sinc(qy*W/2).*sinc(qz*H/2).*exp(-j*qz*H/2);