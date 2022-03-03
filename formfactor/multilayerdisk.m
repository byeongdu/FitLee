function y = multilayerdisk(q, core_radius, core_thick, shell_thick, density, I0, back)
% y = multilayerdisk(q, core_radius, core_thick, shell_thick, density, I0, back)
% ex.
% if there is 3 layer (core, shell1, shell2)
% shell_thick = [shell1, shell2], total 2 members;
% density = [density of core, of shell1, of shell2, of solvent], total 4 member2.

y = zeros(size(q));
V = 1;
dangle = pi/2/64;
angle = 0:dangle:pi/2;
shell_thick = [0, shell_thick];
y1 = y;
yall = y;
ttt = -0.5:0.1:0.5;
nR = core_radius - core_radius*ttt;
rdist = schultzdist(nR, core_radius, 0.1*core_radius);

for k = 1:length(nR)
   for j=1:length(angle)
	alpha = angle(j);
	R = nR(k);
	halfL = core_thick/2;
	y1 = zeros(size(q));
	for i=1:length(shell_thick)
		%R = R + shell_thick(i);
		halfL = halfL + shell_thick(i);
		y1 = y1 + (density(i)-density(i+1))*disk(q, R, halfL, alpha);
	end
	V = vol(R, halfL);
	y = y + I0*sin(alpha)*abs(y1).^2/V*dangle+back;
   end
   yall = yall + rdist(k)*y/length(nR);
end

function y = sisaxs(x)
    y = sinc(x/pi);

function V = vol(R, halfL)
	V = pi*R.^2.*halfL*2;

function y = disk(q, R, halfL, alpha)
    if nargin == 3
	alpha = 0;
    end
    y = 2*vol(R, halfL).*sisaxs(q.*halfL.*cos(alpha)).*besseljc(q.*R.*sin(alpha));
%	y = 2*vol(R, halfL).*sisaxs(q.*halfL.*cos(alpha))./q;
