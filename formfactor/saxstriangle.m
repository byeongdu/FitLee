function F = saxstriangle(qx,qy,R,angle)
%F = saxstriangle(qx,qy,R,angle)
% analytical form factor of a isosceles triangle
% R is the length of the perpendicular line between the acute and the bottom.
% R = A to D.
%            . A
%           .  .
%          .    .
%         .       .
%        .    D    .
%      B. . . . . . .C
%
%
% 'angle' indicate an angle between the two equilaterals in radian.
% the tip is on the coordination center.
% bottom of the triangle is parallel to x axis
% its center of mass is on negative side of y axis.
% derived by Byeongdu Lee. 2006, 2, 12.
j = sqrt(-1);
    if nargin<4
        angle = pi/6;  % equilateral triangle
    end
    u = qx+eps;v=qy+eps;
    L = R*tan(angle);
    %F = -R/u/(u^2*L^2-v^2*R^2)* ...
    %    (-2*u*L + u*L*exp(-i*(u*L-v*R)) + v*R*exp(-i*(u*L-v*R)) + u*L*exp(i*(u*L+v*R)) - v*R*exp(i*(u*L+v*R)));
    U = u*L;V=v*R;expV = exp(j*V);
    F = 2*R./u./(U.^2-V.^2).*(U-U.*cos(U).*expV+V*j.*sin(U).*expV);
