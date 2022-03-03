function y=lorza3D(x,y,z,p)
% lorza3D      : Lorentz Area
% y=lorza2D(x,y,z,p)
% http://en.wikipedia.org/wiki/Lorentzian_function
% p = [ Area x0 y0 z0 Width ]

% Author: B. Lee
% Description:  Lorentzian
% when p(1) =1, the area should be 1.
pcx = p(2);pcy = p(3);pcz = p(4);
c = [pcx, pcy, pcz];
realc = ~isnan(c);
if sum(realc) == 3
    y = p(1)/(2*pi)*p(5)./ (p(5)^2 + (x-pcx).^2 + (y-pcy).^2 + (z-pcz).^2).^1.5;
else 
    q = {x, y, z};
    q = q(realc);
    c = c(realc);
    if sum(realc) == 2
        y = lorza2D(q{1}, q{2}, [p(1), c, p(5)]);
    else
        y = lorza(q{1}, [p(1), c, p(5), 0]);
    end
end