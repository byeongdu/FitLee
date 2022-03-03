function Iq = multilayersphere3(qn, radius, density, sig)
% [y, y1] = multilayersphere3(q, radius, density, sig)
% Here, position of cores are randomly located in the shell..
% Otherwise, identical to multilaersphere2.m
%y = zeros(size(q));
V = 0;
q = linspace(0.0001, 1, 1000)';
%Aq = y;
%Iq = y;
rdiff = radius(2) - radius(1);
corePosition = rdiff*rand(10,1);
%Iq = 0;
Iq = zeros(size(q));
for j=1:numel(corePosition)
    Aq = zeros(size(q));
    for i=1:length(radius)
        Aq = Aq + (density(i)-density(i+1))*sphereamp(q, radius(i));
        if i==1
            Aq = Aq.*exp(sqrt(-1)*q*corePosition(j));
        end
        V = V + (density(i)-density(i+1))*radius(i).^3/3*4*pi;
    end
yt = abs(Aq/V).^2;
Iq = Iq + yt*V^2;
end
Iq = formfactorscale(q, Iq, radius(end), qn, radius(end), sig);

