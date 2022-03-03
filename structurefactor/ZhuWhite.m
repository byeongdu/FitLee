function [y, name, pnames, pin] = ZhuWhite(x, p, flag)
% [y, name, pnames, pin] = ZhuWhite(x, p, flag)
% Zhu et al. J. Chem. Phys., 1995, 104(22), 1996
% p =('a0','b', 'Df' 'r0', 'Contrast^2', 'bkg')
if nargin == 2;
    a0  = p(1);
    b = p(2);
    Df = p(3);
    K = p(4);
    bkg = p(5);
    a0 = abs(real(a0)); b= abs(real(b)); Df = abs(real(Df)); K = abs(K);
    sinterm1 = sin((Df-1)*atan((x+b)*a0)).*(1+(x+b).^2*a0^2).^(1/2-Df/2);
    sinterm2 = sin((Df-1)*atan((x-b)*a0)).*(1+(x-b).^2*a0^2).^(1/2-Df/2);
    y = 2*pi*a0^Df*K./(a0*x).*(sinterm1 + sinterm2).*gamma(Df-1);
    y = y+bkg;

else
	y=[];
	name='Gauss distr Sphere';
	pnames=str2mat('a0','b', 'Df', 'Contrast^2', 'bkg');
	if flag==1, pin=[74 0.024 2.6 1 0]; else pin = p; end
end
