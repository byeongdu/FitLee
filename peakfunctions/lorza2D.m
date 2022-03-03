function y=lorza2D(x,y,p)
% lorza2D      : Lorentz Area
% y=lorza2D(x,y,p)
% http://en.wikipedia.org/wiki/Lorentzian_function
% p = [ Area Centrex centery Width ]

% Author: B. Lee
% Description:  Lorentzian
% when p(1) =1, the area should be 1.


c = p(2:3);
realc = ~isnan(c);
if sum(realc) == 2
    y = p(1)/(2*pi)*p(4)./ (p(4)^2 + (abs(x)-p(2)).^2 + (y-p(3)).^2 ).^1.5;
else 
    q = {abs(x), y};
    q = q(realc);
    c = c(realc);
    y = lorza(q{1}, [p(1), c, p(4), 0]);
end