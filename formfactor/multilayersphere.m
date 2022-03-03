function [Iq, Aq, y] = multilayersphere(q, radius, density, sig)
% [y, y1] = multilayersphere(q, radius, density, sig)
V = 0;
if (nargin == 3)
	sig = 0;
end
qn = q;
q = linspace(0.0001, 1, 1000)';
y = zeros(size(q));

%if (sig ~= 0)
%	%r = 1:radius(1)/20:(radius(1)*3);
%    r = linspace(1, radius(1)*3, 20);
%	%nr = logn(r, [radius(1), sig, 1]);
%    nr = schultzdist(r, radius(1), sig);
%else
	nr = 1;
	r = radius(1);
%end

%        M2 = 0;
    
%        for i = 2:NumS
%            M2 = M2 + (4/3*pi*Radi(i)^3)*(rho(i)-rho(i-1));
%        end
        
%        M = rho(1)*(4/3*pi*Radi(1)^3) + M2;
Aq = y;
Iq = y;
for j=1:length(nr)
	rad = radius/radius(1)*r(j);
    V = 0;
    yt = 0;
    Aq = 0;
	for i=1:length(rad)
		Aq = Aq + (density(i)-density(i+1))*sphereamp(q, rad(i));
		V = V + (density(i)-density(i+1))*rad(i).^3/3*4*pi;
    end
    yt = nr(j)*abs(Aq/V).^2;
    Iq = Iq + yt*V^2;
	y = y + yt;
    %y = y + nr(j)*abs(y1).^2;
end
y = y/sum(nr);
Iq = Iq/sum(nr);
Iq = formfactorscale(q, Iq, radius(end), qn, radius(end), sig);
if numel(nr) > 1
    Aq = [];
end