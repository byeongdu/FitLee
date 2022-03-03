function y = background2d(x, y, p)
% 2D background with 2nd polynomial.
% bkg = p(1) + p(2)*x + p(3)*y;
% bkg = p(1) + p(2)*x + p(3)*y + p(4)*x.^2 + p(5)*y.^2;

if numel(p) == 5
    c = p(1);
    c1 = p(2);
    c2 = p(3);
    c3 = p(4);
    c4 = p(5);
    y = c+c1*x+c2*y+c3*x.^2+c4*y.^2;
end
if numel(p) == 3
    c = p(1);
    c1 = p(2);
    c2 = p(3);
    y = c+c1*x+c2*y;
end