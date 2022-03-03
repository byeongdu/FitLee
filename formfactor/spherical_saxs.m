function [y, name, pnames, pin] = spherical_saxs(q, p, flag)
% [y, name, pnames, pin] = spherical_saxs(x, p, flag)
global xx;
global eden;

if nargin == 2;
    amp = p(1:2:11);
    sig = abs(p(2:2:12));
    dmax = p(13);
    I0 = p(14);
    back = p(15);
    r = 0:dmax/6:dmax;
    r(end) = [];
    [xx, eden] = eden_line(r, amp, sig, dmax);
    eden(end) = 0;
    eden = [eden, 0];
    y = multilayersphere(q, xx, eden);
    y = I0*abs(y).^2 + back;
    sig(2:end) = vect2row(sig(2:end))./vect2row(r(2:end));
%    y = y.*q.^4;
%    kk = smearwl([q, y], [12, 0.6]);y = kk(:,2).*q.^2;
    pin(1:2:11) = amp;
    pin(2:2:12) = sig;
else
	y=[];
	name='multilayered spherical shell';
	pnames=str2mat('amp1', 'sig1', 'amp2', 'sig2', 'amp3', 'sig3', 'amp4', 'sig4', 'amp5', 'sig5', 'amp6', 'sig6', 'dmax', 'I0', 'back');
	if flag==1, pin=[-1, 0.1, -1, 0.1, 0, 0.1, 1, 0.1, 2, 0.1, 0.5, 0.1, 50, 1, 0]; else pin = p; end
end

function [xx, eden] = eden_line(r, amp, sig, dmax)
sig1 = sig(1);
    r = vect2row(r);
    sig = vect2row(sig);
    t = find(sig>1);
%    sig(t) = 0.3;
sig(1) = sig1;
    sig(2:end) = sig(2:end).*r(2:end);
    xx = 0:dmax/25:dmax;
    eden = zeros(size(xx));
    amp = amp./abs(amp(end));
    for i=1:length(r)
        [xt, yt] = gaus(xx, r(i), sig(i));
        eden = eden + amp(i)*yt';
    end
