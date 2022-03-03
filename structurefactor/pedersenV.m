function [y, name, pnames, pin] = pedersenV(x, p, flag)
% [y, name, pnames, pin] = pedersen(x, p, flag)
% Pedersen formalism to calculate scattering intensity.
% Local monodisperse approximation assumed.
% Formfactor assumed as sphere scattering.
% size distribution function : Gaussian(0), Gamma(1)
% Structure factor : Hard sphere
if nargin == 2;
    center  = p(1);
    width = p(2);
    I0 = p(3);
    bkg = p(4);
    volfrac = p(5);
    distrFunc = p(6);
    dmax = p(7);
    
    if distrFunc == 0
        if dmax == 0
            [gDist, xx] = gaussdist99(center, width);
        else
            [gDist, xx] = gaussdist99(center, width, dmax);
        end
    elseif distrFunc == 1
        [gDist, xx] = gamdist99(center, width);
    elseif distrFunc == 2
        if dmax == 0
            [gDist, xx] = logndist99(center, width);
        else
            [gDist, xx] = logndist99(center, width, dmax);
        end
    end

    [x1, x2]= size(x);
    
    if x1 < x2
        x = x';
    end
    %normFactor = sum(sum((gDist*xx.^6)));size(normFactor);%SS = SFactor(x, xx, volfrac);
    %y = sum(spheretype(x, xx).*SS*(gDist*xx.^6), 2)*I0/normFactor + bkg ;
    normFactor = sum(gDist);
    y = zeros(size(x));
%    figure;hold on
    for i = 1:length(xx)
        Form = spheretype(x, xx(i));
        Struc = SFactor(x, xx(i), volfrac);
        y = y + Form.*Struc.*gDist(i);
%        plot(x, Form, 'b');plot(x, Struc, 'r');
    end
    y = y*I0/normFactor+bkg;

else
	y=[];
	name='Gauss distr Sphere';
	pnames=str2mat('center','width', 'I0', 'bkg', 'Volume fraction', 'Distribution Function','dmax');
	if flag==1, pin=[10 1 1 0, 0.1 0, 0]; else pin = p; end
end

function SF = SFactor(q, Rhs, vf)
% vf is volume fraction of hard sphere
alpha = (1 + 2*vf)^2/(1-vf)^4;
beta = -6*vf*(1 + vf/2)^2/(1-vf)^4;
gamma = vf*alpha/2;

if length(Rhs) > 1
    [R, q] = meshgrid(Rhs, q);
    A = R.*q*2;
else
    A = Rhs.*q*2;  
end

GA = alpha*(sin(A)-A.*cos(A))./A.^2 + ...
    beta*(2.*A.*sin(A) + (2 - A.^2).*cos(A)- 2)./A.^3 + ...
    gamma*(-1*A.^4.*cos(A) + 4*((3*A.^2-6).*cos(A) + (A.^3 - 6 * A).*sin(A) + 6))./A.^5;
SF = 1./(1+24*vf*GA./A);
%plot(SF(:,1))