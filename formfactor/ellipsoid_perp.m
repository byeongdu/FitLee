function Formfactor = ellipsoid_perp(qy, qz, H, W, xyangle)
% an oblate is perpendicular to x-ray beam (x-axis).
% short axis (H) is x-axis
% two long axes (W's) are y and z axis...

R = W/2;

        % rotate pyramid 45degree in xy plane.
    
np = 100;
dx = H/np;
xdir = linspace(0, H/2, np);
%lengthx = np+1;


y = qy;
x = zeros(size(qy));

xyangle = deg2rad(xyangle);
qx = cos(xyangle)*x + sin(xyangle)*y;
qy = -sin(xyangle)*x + cos(xyangle)*y;

%x = abs(xdir - H + R);
%Rx = sqrt(R.^2 - x.^2);
    

F = zeros(size(qy));
for i=1:np;
    x = xdir(i);
    Rx = R*sqrt(1-4*x^2/H^2);
    qp = sqrt(qy.^2+qz.^2);
    qpRx = qp.*Rx;
    [xzero, yzero] = find(qpRx == 0);
    qpRx(xzero, yzero) = 0.00000001;
    Fb = 4*pi*Rx.^2.*bessel(1, qpRx)./(qpRx).*cos(qx*x)*dx;
    F = F+Fb;
end
Formfactor = F;
%Formfactor = sum(F, 3);
