function y = sinc(x)
% sinc function
% sinc = sin(x)/x
% Byeongdu Lee
t = find(abs(x) > eps);
y = ones(size(x));
y(t) = sin(x(t))./x(t);
end