function [y, name, pnames, pin]=cubetype(x, p, flag)
% Form factor of Cube or rectangular parallelpipedons
% variable, x, p, flag
% p = [a b c]

a = p(1);b=p(2);c=p(3);

AL = pi/50:pi/50:pi/2;
BE = pi/50:pi/50:pi/2;
[q, al, be] = meshgrid(x, AL, BE);

y = sin(q.*a.*sin(al).*cos(be))./(q.*a.*sin(al).*cos(be)).*sin(q.*b.*sin(al).*cos(be))./(q.*b.*sin(al).*sin(be)).*sin(q.*c.*cos(al))./(q.*c.*cos(al)).*sin(al);
y = sum(sum(y, 3), 1);
y = y/y(1);