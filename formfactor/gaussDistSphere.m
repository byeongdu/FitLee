function [y, name, pnames, pin] = gaussDistSphere(x, p, flag)

if nargin == 2;
    center  = p(1);
    width = p(2);
    I0 = p(3);
    bkg = p(4);

    [gDist, xx] = gaussdist99(center, width);

    [x1, x2]= size(x);
    
    if x1 < x2
        x = x';
    end
    normFactor = sum(sum(gDist*xx'.^6));
    intensity = spheretype(x, xx)*(gDist*xx'.^6);%figure(100);loglog(x, spheretype(x, xx))
    y = sum(intensity, 2)*I0/normFactor + bkg ;

else
	y=[];
	name='Gauss distr Sphere';
	pnames=str2mat('center','width', 'I0', 'bkg');
	if flag==1, pin=[10 1 1 0]; else pin = p; end
    
end