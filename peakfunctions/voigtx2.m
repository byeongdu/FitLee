function [y, name, pnames, pin]=gaussx2(x,p, flag)
% voigtx2   : 2 voigts
% [y, {name, pnames, pin}]=voigtx2(x,p, {flag})
%
% MFIT 2 Voigt fitting functions
% p = [ Amplitude Centre Gauss_width Lorz_width Amplitude2 Centre2 Gauss_width2 Lorz_width2 Background ]

% Author:  BLee <leedo@postech.ac.kr>
% Description:  2 voigts

global ngss;

if nargin==2
    N = 32;
	b = -sqrt(log(2))/p(3);
	a = b*p(4);
	b = b*2*i;
	z = a + b*(x-p(2));
    
    b2 = -sqrt(log(2))/p(7);
	a2 = b2*p(8);
	b2 = b2*2*i;
	z2 = a2 + b2*(x-p(6));

	M=2*N; M2=2*M; k=[-M+1:1:M-1]';
	L=sqrt(N/sqrt(2));
	tt=(L*tan(k*pi/M2)).^2;
	f=[0; exp(-tt).*(L^2+tt)];
	a=real(fft(fftshift(f)))/M2;
	a=flipud(a(2:N+1));

    l=L-z;
	Z=(L+z)./l;
	pp=polyval(a,Z);
	y1=p(1)*real(2*pp ./l.^2+(1/sqrt(pi))*ones(size(z)) ./l);

    l2=L-z2;
	Z2=(L+z2)./l2;
	pp2=polyval(a,Z2);
	y2=p(5)*real(2*pp2 ./l2.^2+(1/sqrt(pi))*ones(size(z2)) ./l2);

    y = y1+y2 + p(9)*x+p(10);
else
	y=[];
	if flag==1, pin=[0 0 0 1 0 0 0 1 0 0]; else pin = p; end
	name='2 Gaussians';
	pnames=str2mat('Amplitude1', 'Centre1', 'Gauss_width1', 'Lorz_width1',...
                   'Amplitude2', 'Centre2', 'Gauss_width2', 'Lorz_width2', 'Background1', 'back2');

   if flag==2

      mf_msg('Click on background');    
		[x bg]=ginput(1);

        mf_msg('Click on peak1');
		[cen1 amp]=ginput(1);
		mf_msg('Click on width1');
		[width1 y]=ginput(1);
		width1=abs(width1-cen1);
		amp1=amp-bg;

        mf_msg('Click on peak2');
		[cen2 amp]=ginput(1);
		mf_msg('Click on width2');
		[width2 y]=ginput(1);
		width2=abs(width2-cen2);
		amp2=amp-bg;

        pin=[2*amp1 cen1 width1/2 width1/2  2*amp2 cen2 width2/2 width2/2  bg, bg2];
	end

end
