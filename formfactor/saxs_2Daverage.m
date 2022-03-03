function y = saxs_2Daverage(q, funcname, parameters)
strcmd = ['F = @', funcname, ';'];
eval(strcmd);
dtheta = pi/100;
%dphi = pi/100;
theta = 0:dtheta:pi/2;
%phi = 0:dphi:pi/2;
[q, theta] = ndgrid(q, theta);
qx = q.*cos(theta);
qy = q.*sin(theta);
qz = 0;
y = abs(feval(F, qx, qy, qz, parameters)).^2;
y = trapz(y, 2)*dphi*2/pi;