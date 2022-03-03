function [y, name, pnames, pin, gr] = ZhuStruc(x, p, flag)
% [y, name, pnames, pin] = ZhuStruc(x, p, flag)
% Zhu et al. J. Chem. Phys., 1995, 104(22), 1996
% pnames=str2mat('a0','b', 'Df', 'Contrast^2', 'bkg');
if nargin == 2;
    w  = p(1);
    D = p(2);
    K = p(3);
%    r0 = p(3);
    r = eps:1:50; r=r/50*0.5*D;lenr = length(r);
%    k = find(r<r0);r(k) = eps;
    [xx, yx] = size(x);
    if xx<yx
        x = x';
    end
    
    dr = 0.5/200*D;
    [q2, r2] = meshgrid(x, r);
    gr = 1-exp(-r2/w).*sin(2*pi*r2/D)./(2*pi*r2/D);
    %gr = -exp(-r2/w).*cos(2*pi*r2/D);
    rq = besselj(0, q2.*r2);
    grr = (gr-1).*rq.*r2*dr;
    y = K*pi*sum(grr, 1)/lenr/max(r)'+1;
else
	y=[];
	name='Zhu Structure factor';
	pnames=str2mat('w','P', 'K');
	if flag==1, pin=[10 100, 15]; else pin = p; end
end
