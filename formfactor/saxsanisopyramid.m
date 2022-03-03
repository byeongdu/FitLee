function [y, V] = saxsanisopyramid(qx, qy, qz, parameter)
% referece: G. Renaud et al. Surface Science Reports 64(2009) 255-380
% Bottom plane of the anisotropic pyramid is a retangle, whose
% 2R is the short length and along the x axis.
% 2W is the long axis and along the y axis.
% H is the height along the z axis.
% The center of a mass of the bottom retangle is at (0,0,0)
%


if numel(parameter) >= 3
    R = parameter(1);
    W = parameter(2);
    H = parameter(3);
else
    error('the number of parameters should be at least 3 for saxsanisopyramid.m')
end
if numel(parameter) == 3
    alpha = atan(H/R);
end
i = sqrt(-1);
V = 4*(W*R*H-H^2*(R+W)/(2*tan(alpha))+H^3/(3*tan(alpha)^2));
qx = -qx+eps;
qy = -qy+eps;
qz = -qz;
t = find((abs(qx) < 1000*eps) & (abs(qy) < 1000*eps));
qx(t) = 1E-7;
qy(t) = 1E-8;
q1 = 1/2*((qx-qy)/tan(alpha)+qz);
q2 = 1/2*((qx-qy)/tan(alpha)-qz);
q3 = 1/2*((qx+qy)/tan(alpha)+qz);
q4 = 1/2*((qx+qy)/tan(alpha)-qz);
K1 = sinc(q1*H).*exp(i*q1*H)+sinc(q2*H).*exp(-i*q2*H);
K2 = -i*sinc(q1*H).*exp(i*q1*H)+i*sinc(q2*H).*exp(-i*q2*H);
K3 = sinc(q3*H).*exp(i*q3*H)+sinc(q4*H).*exp(-i*q4*H);
K4 = -i*sinc(q3*H).*exp(i*q3*H)+i*sinc(q4*H).*exp(-i*q4*H);

y = H./(qx.*qy).*(cos(qx*R-qy*W).*K1 ...
    +sin(qx*R-qy*W).*K2 ...
    -cos(qx*R+qy*W).*K3 ...
    -sin(qx*R+qy*W).*K4);
