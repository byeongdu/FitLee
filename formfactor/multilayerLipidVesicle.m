function y = multilayerLipidVesicle(q, radius, density, sig, I0, back)
y = zeros(size(q));
V = 0;
if (nargin == 3)
	sig = 0;
	I0 = 1;
	back = 0;
end

if (nargin == 4)
	I0 =1;
	back = 0;
end

if (sig ~= 0)
	r = (radius(1)-radius(1)*0.8):radius(1)/20:(radius(1)+radius(1)*0.8);
	%nr = schultzdist(r, radius(1), sig, 1);
	nr = logn(r, [radius(1), sig, 1]);
else
	nr = 1;
	r = radius(1);
end

y1 = y;
for j=1:length(nr)
	rad = radius - radius(1)+r(j);
	for i=1:length(rad)
		y1 = y1 + (density(i)-density(i+1))*sphereAmp(q, rad(i));
		V = V + (density(i)-density(i+1))*rad(i).^3/3*4*pi;
	end
	y = y + I0*nr(j)*abs(y1).^2+back;
end
