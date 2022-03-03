function [y, name, pnames, pin] = AsymDGaussCum(x, p, flag)
% function [y, name, pnames, pin] = AsymDGaussCum(x, p, flag)
% Asymetric Double Gauss Cumulative function.

if nargin == 2;

    y = p(1)/2*(1+erf((x-p(2)+p(3)/2)/(sqrt(2)*p(4)))).*(1/2-1/2*erf((x-p(2)-p(3)/2)/(sqrt(2)*p(5))));
    
    else
	y=[];
	name='Gauss distr Sphere';
   
    pnames = str2mat('MaxAmp','Center', 'Width', 'shape1', 'shape2');
	if flag==1, pin=[1, 0, 0.1, 0.1, 0.1]; else pin = p; end
end
