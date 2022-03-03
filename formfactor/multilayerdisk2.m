function yall = multilayerdisk2(q, core_radius, layerposition, density, sigofR, I0, back)
% y = multilayerdisk(q, core_radius, core_thick, shell_thick, density, I0, back)
% ex.
% if there is 3 layer (core, shell1, shell2)
% shell_thick = [shell1, shell2], total 2 members;
% density = [density of core, of shell1, of shell2, of solvent], total 4 member2.

y = zeros(size(q));
V = 1;
dangle = pi/2/64;
angle = 0:dangle:pi/2;
y1 = y;
yall = y;

%angleinput = input('Enter angle : ');
%if ~isempty(angleinput)
%    angle = angleinput;
%end

if (sigofR > 0)
    ttt = -0.5:0.1:0.5;
    nR = core_radius - core_radius*ttt;
    rdist = schultzdist(nR, core_radius, sigofR);
else
    nR = core_radius;
    rdist = 1;
end

layerthick = layerposition(2)-layerposition(1);
angle = pi/2;

for k = 1:length(nR)
   for j=1:length(angle)
	alpha = angle(j);
	R = nR(k);
	y1 = zeros(size(q));
    
    y1 = y1 + density(1).*sisaxs(layerthick/2*q);
	for i=2:length(density)
		y1 = y1 + 2*density(i).*sisaxs(layerthick/2*q).*cos(layerposition(i)*q);
    end
    
	y = y + I0*sin(alpha)*abs(y1).^2/V*dangle+back;
   end
   yall = yall + rdist(k)*y/length(nR);
end

function y = sisaxs(x)
if (x~=0)
    y = sin(x)./x;
else
    y = 1;
end
%    y = sinc(x/pi);

function V = vol(R, halfL)
	V = pi*R.^2.*halfL*2;

function y = disk(q, R, halfL, alpha)
    if nargin == 3
	alpha = 0;
    end
    y = 2*vol(R, halfL).*sisaxs(q.*halfL.*sin(alpha)).*besseljc(q.*R.*cos(alpha));
%	y = 2*vol(R, halfL).*sisaxs(q.*halfL.*cos(alpha))./q;
