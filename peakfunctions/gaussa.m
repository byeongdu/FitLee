function [y, name, pnames, pin]=gaussa(x,p, flag)
% gauss     : Gaussian area
% [y, {name, pnames, pin}]=gaussa(x,p, {flag}) 
% MFIT Gaussian fitting function
% p = [ Area Centre Width BackG ]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description: Gaussian

if nargin==2;
    y=p(4)+p(1)/sqrt(2*pi)/p(3)*exp(-0.5*((x-p(2))/p(3)).^2);
else
	y=[];
	name='gaussian area';
	pnames=str2mat('area','Centre','Width','Background');
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
