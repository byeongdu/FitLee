function [y, name, pnames, pin]=lorza(x,p, flag)
% lorz      : Lorentz Area
% [y, {name, pnames, pin}]=lorz(x,p, {flag})
%
% MFIT Lorentzian fitting function
% p = [ Area Centre Width Background ]
% integral is : pi*p(1)*p(3) when bg = 0

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description:  Lorentzian

% http://mathworld.wolfram.com/LorentzianFunction.html
% FWHM = p(3)*2
if nargin==2
    FWHM = p(3);
    halfgam = FWHM/2;
	y=p(4)+p(1)/pi/halfgam./ (1 + ((x-p(2))/halfgam).^2 );
else
	y=[];
	name='Lorentz Area';
	pnames=str2mat('Area','Centre','Width','Background');
	if flag==1, pin=[0 0 1 1]; else pin = p; end
	if flag==2
		mf_msg('Click on peak');
		[cen amp]=ginput(1);
		mf_msg('Click on width');
		[width y]=ginput(1);
		width=abs(width-cen);
		mf_msg('Click on background');
		[x bg]=ginput(1);
		amp=amp-bg;
		pin=[amp cen width bg];
	end
end
