function [y, name, pnames, pin] = pedersen(x, p, flag)
% [y, name, pnames, pin] = pedersen(x, p, flag)
% Pedersen formalism to calculate scattering intensity.
% Local monodisperse approximation assumed.
% Formfactor assumed as sphere scattering.
% size distribution function : Gaussian(0), Gamma(1), lognormal(2),
% schultz(3)
% Structure factor : Hard sphere
%    center  = p(1);
%    width = p(2);
%    I0 = p(3);
%    bkg = p(4);
%    volfrac = p(5);
%    distrFunc = p(6);
%    dmax = p(7);
k = find(x == 0);
if ~isempty(k)
    if (k>0) & (k<length(x))
        x(k) = x(k+1);
    end
end


if nargin == 2;
    center  = p(1);
    width = p(2);
    I0 = p(3);
    bkg = p(4);
    volfrac = p(5);
    distrFunc = p(6);
    dmax = p(7);
%    power = p(10);
%    pf = p(11);
%    if length(p) > 9
%        pos = p(8);
%        pw = p(9);
%    else
%        pos = 0;
%        pw = 0;
%    end
    
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
            [gDist, xx] = logndist9(center, width);
        else
            [gDist, xx] = logndist9(center, width);
        end
    end

    [x1, x2]= size(x);
    
    if x1 < x2
        x = x';
    end
    %normFactor = sum(sum((gDist*xx.^6)));size(normFactor);%SS = SFactor(x, xx, volfrac);
    %y = sum(spheretype(x, xx).*SS*(gDist*xx.^6), 2)*I0/normFactor + bkg ;
    if distrFunc <3
        normFactor = vect2row(gDist)*vect2row(xx.^6)';
    elseif distrFunc == 3
        normFactor =1;gDist=1;
        xx = center;
        Form = SchultzsphereFun(x, center, width);
        Struc = SFactor(x, center*dmax, volfrac);
    end
    y = zeros(size(x));
%    figure;hold on
    for i = 1:length(xx)
        if distrFunc < 3
                Form = spheretype(x, xx(i))*xx(i)^6;
            if (volfrac == 1)
                Struc = SFactor2Dpara(x, pos, pw);
            else
                Struc = SFactor(x, xx(i), volfrac);
            end
        end

        y = y + Form.*Struc.*gDist(i);
%        plot(x, Form, 'b');plot(x, Struc, 'r');
    end
    y = y*I0/normFactor+bkg;% + pf*x.^power;

else
	y=[];
	name='Gauss distr Sphere';
	pnames=str2mat('center','width', 'I0', 'bkg', 'Volume.fraction', 'Distribution.Function','dmax', 'position', 'positionw', 'power', 'prefactor');
	if flag==1, pin=[10 1 1 0, 0.1 0, 0, 20, 0.1, -2, 1]; else pin = p; end
end

function SF = SFactor(q, Rhs, vf)
% vf is volume fraction of hard sphere
if vf == 0 
    SF = ones(size(q));
    return
end

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

function SF = SFactor2Dpara(q, dsp, sig)
% vf is volume fraction of hard sphere
F = exp(-1*q.^2*sig^2);
SF = (1-F.^2)./(1-2*F.*cos(q*dsp)+F.^2);
%plot(SF(:,1))
