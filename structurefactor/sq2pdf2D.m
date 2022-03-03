function y = sq2pdf2D(q, sq, r, n_p)
% function y = sq2pdf(q, sq, r, n_p)
% for 2D.
% np*(gr-1) = 1/pi^2*int(q*iq*sin(qr))dq
% iq = sq-1 for a single particle system
% or sq = 1+beta*iq for multi-particle system.

t = isinf(sq);
sq(t) = [];q(t)=[];
dq = q(2)-q(1);
%dsp = 0.2732;
%R = dsp/2;

y = zeros(size(r));
for i=1:length(r)
%    y(i) = 1/pi^2*sum(q.*(sq-1).*sin(q*r(i)))*dq/n_p+1;  % for single componentsystem
    y(i) = 1/(2*pi)/n_p*sum(1./r.*(sq-1).*sin(q*r(i)))*dq+1;  % for single componentsystem
end
