function y = besseljc(x)
y = besselj(1, x+eps)./(x+eps);