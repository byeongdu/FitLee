function y=saxsparallelpiped(qx, qy, qz, xl, yl, zl, flag)
% Form factor of Cube or rectangular parallelpipedons
% whose center of mass is at [0,0,0]
% variable, x, p, flag
% p = [a b c]
y = sinc(xl*qx/pi).*sinc(yl*qy/pi).*sinc(zl*qz/pi);
