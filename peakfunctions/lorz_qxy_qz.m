function I=lorz_qxy_qz(x,y,p)
% lorentz function along the q direction.
% gaussian function for the orientational distribution.
% y=lorz_q(x,y,p)
% http://en.wikipedia.org/wiki/Lorentzian_function
% p = [ Area Centrex centery width_qz width_qxy ]

% Author: B. Lee
% when p(1) =1, the area should be 1.

I1 = lorza(abs(x), [p(1), p(2), p(5), 0]);
Ig = lorza2D(abs(x), y, [1, p(2), p(3), p(4), 0]);
% calculate 2D diffraction with broader width, which is width_qz
% and calculate 1D diffraction with narrower width, which is width_qxy
% multiply narrow one to the 2D.
% so, along the qxy, it becomes narrower.

I = I1.*Ig;