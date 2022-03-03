function y=SchultzCoreShellFF2(q, p)
% Analytical solution for a polydisperse spherical core-shell particles.
% Only core is polydisperse
% reference : P. Barlett and R. H. Ottewill, J. Phy. Chem. (1992). 96, 3306-3318

r = p(1);	   % average radius Rshell
pd = p(2);	   % polydisperse pd = 1/sqrt(z+1), which is difference from standard deviation.
t = p(3);      % thickness of shell
rho1 = p(4);   % rho1 = rho_core
rho2 = p(5);   % rho2 = rho_shell
rho3 = p(6);   % rho3 = rho_solv

if pd == 0
   if (beta == 1) | (beta == 0)
    	y = spheretype(q, r);
   else
	y = SphConcShell(q, [2, rho1, beta*r, rho2, r, 1]);
   end
else
%	z=r*r/(fwhm*fwhm);
	z = 1/pd^2-1;
	x = q*r;
	y = q*t;
	gamma = (rho3-rho2)/(rho1-rho2);
	Bx = Bx(z, x);
	Dx = Dx(z, x);
	c1= c1(y, gamma);
	c2= c2(y, gamma);
	c3= c3(y, gamma);
	c4= c4(y, gamma, c1);
	c5= c5(y, gamma, c2);
	c6= c6(y, gamma, c3);
	c7= c7(y, gamma, c5);
	c8= c8(y, gamma, c4);
	c9= c9(y, gamma);
	a = 16*pi^2./q.^6.*(rho2-rho1)^2;
	b = (c1 + c2.*x + c3.*x.^2.*((z+2)./(z+1))...
		+ Bx.^((z+1)/2).*(c4.*cos((z+1).*Dx) + c7.*sin((z+1).*Dx))...
		+ x.*Bx.^((z+2)/2).*(c5.*cos((z+2).*Dx) + c8.*sin((z+2).*Dx))...
		+ (z+2)./(z+1).*x.^2.*Bx.^((z+3)/2).*(c6.*cos((z+3).*Dx) + c9.*sin((z+3).*Dx)));
	y = a.*b;
end

function y = c1(y, gamma)
	y= 1/2-gamma*(cos(gamma)+y.*sin(y)) + gamma^2/2*(1+y.^2);

function y = c2(y, gamma)
	y= gamma*y.*(gamma-cos(y));

function y = c3(y, gamma)
	y= (gamma^2+1)/2-gamma*cos(y);

function y = c4(y, gamma, c1)
	y= gamma^2*(y.*cos(y)-sin(y)).^2-c1;

function y = c5(y, gamma, c2)
	y= 2*gamma*sin(y).*(1-gamma*(y.*sin(y)+cos(y)))+c2;

function y = c6(y, gamma, c3)
	y= c3-gamma^2*sin(y).^2;

function y = c7(y, gamma, c5)
	y= gamma*sin(y)-gamma^2/2*(1+y.^2).*sin(2*y)-c5;

function y = c8(y, gamma, c4)
	y= c4-1/2+gamma*cos(y)-gamma^2/2*(1+y.^2).*cos(2*y);

function y = c9(y, gamma)
	y= gamma*sin(y).*(1-gamma*cos(y));

function y = Bx(z, x)
	y = (z+1)^2./((z+1)^2+4*x.^2);

function y = Dx(z, x)
	y = atan(2*x/(z+1));
