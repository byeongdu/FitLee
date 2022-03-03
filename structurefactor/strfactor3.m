function [y, name, pnames, pin] = strfactor3(q, p, flag)
% function [y, name, pnames, pin] = strfactor3(q, p, flag)
% p = [dsp, sig];
% 2D paracrystal
if nargin==2;
    dsp = p(1);
    sig = p(2);
    F = exp(-1*q.^2*sig^2);
    y = (1-F.^2)./(1-2*F.*cos(q*dsp)+F.^2);
    
else
	y=[];
	name='Gaussian_2Dparacrystal';
	pnames=str2mat('position', 'gau_Width_position');
	if flag==1, pin=[20, 0.1]; else pin = p; end
end

%plot(SF(:,1))