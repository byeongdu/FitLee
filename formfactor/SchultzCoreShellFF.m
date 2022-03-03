function [y, C] = SchultzCoreShellFF(q, varargin)
% Analytical solution for a polydisperse spherical core-shell particles.
% fixed core/shell ratio
% reference : Joachim Wagner, J. Appl. Cystal. (2004). 37, 750-756
% y=SchultzCoreShellFF(q, p)
%   p = [R_overall, polydispersity(<1), Rcore/Rshell, Rho_core-rho_solvent,
%   Rho_shell - rho_solvent];
%   or
%   p = [R(array; Rcore, R_overall), edensity(array; core, shell, solvent), polydispersity(<1), ];
%
if numel(varargin) == 5
    r = varargin{1};	   % average radius Rshell
    pd = varargin{2};	   % polydisperse pd = 1/sqrt(z+1), which is different from standard deviation.
    beta = varargin{3};   % beta = Rcore/Rshell; beta = 1 for a homogeneous sphere
    rho1 = varargin{4};   % rho1 = rho_core - rho_solv
    rho2 = varargin{5};   % rho2 = rho_shell - rho_solv
else
    r = varargin{1};
    eden = varargin{2};
    pd = varargin{3};
    beta = r(1)/r(2);
    r = r(2);
    rho1 = eden(1)-eden(3);
    rho2 = eden(2)-eden(3);
end

N = 1000;

if pd == 0
   if (beta == 1) | (beta == 0)
    	y = spheretype(q, r);
   else
	y = SphConcShell(q, [2, rho1, beta*r, rho2, r, 1]);
   end
else
%	z=r*r/(fwhm*fwhm);
	z = 1/pd^2-1;
    if z<= 50
        N = 2;
    elseif z<=120
        N = 10;
    elseif z<=160
        N = 100;
    elseif z>160
        error('polydispersity should be larger, otherwise I cannot calculate')
    end
	GM = (z+1)./(q*r);
	C = Contrast(beta, rho1, rho2);
    
    y2 = (GM).^(z*1/N);
	
    F1p = F1pf(z, beta, GM);
	F1m = F1mf(z, beta, GM);
	F2p = F2pf(z, beta, GM);
	F2m = F2mf(z, beta, GM);
	F3p = F3pf(z, beta, GM);
	F3m = F3mf(z, beta, GM);
	G1 = G1f(z, beta, GM);
	G2 = G2f(z, beta, GM);
	G3 = G3f(z, beta, GM);
	L1 = L1f(z, beta, GM);
	L2 = L2f(z, beta, GM);
	L3 = L3f(z, beta, GM);
	K1 = K1f(z, beta, GM);
	K3 = K3f(z, beta, GM);
    % z! = gamma(z+1)
	%y = 9*(GM).^(z+7)*1/C^2*gamma(z+1)/gamma(z+7));
	y = (rho2*(rho2-rho1)*(F1p-F1m+(1+beta)*F2p+(1-beta)*F2m-beta*(F3p+F3m))...
	    -(rho2-rho1).^2.*(G1/2+beta*G2-beta^2*G3/2) ...
	    - rho2^2*(L1/2+L2-L3/2) ...
	    + (rho1^2/2+rho2^2-rho1*rho2)*K1 ...
	    + (beta^2*(rho2-rho1)^2 + rho2^2)*K3/2);
	y = y.*(9*(GM).^(7)*1/C^2*gamma(z+1)/gamma(z+7));
	%y = y.*(GM).^(z/5);
	%y = y.*(GM).^(z/5);
	%y = y.*(GM).^(z/5);
	%y = y.*(GM).^(z/5);
	%y = y.*(GM).^(z/5);
end

function y = Contrast(beta, rho1, rho2)
	y = 4*pi/3*(rho1.*beta.^3+rho2.*(1-beta.^3));
end
function y = F1pf(z, beta, GM)
	y = cos((z+1).*atan((beta+1)./GM));
    y1 = (GM.^2+(beta+1).^2).^(-(z+1)/2/N);
	%y2 = y.*(GM).^(z*1/N);
    for i=1:N
        y = y.*y1.*y2;
    end
end
function y = F1mf(z, beta, GM)
	y= cos((z+1).*atan((beta-1)./GM));
    y1= (GM.^2+(beta-1).^2).^(-(z+1)/2/N);
	%y2= (GM).^(z*1/N);
    for i=1:N
        y = y.*y1.*y2;
    end
end
function y = F2pf(z, beta, GM)
	y= sin((z+2).*atan((beta+1)./GM)).*(z+1);
    y1= (GM.^2+(beta+1).^2).^(-(z+2)/2/N);
	%y2= (GM).^(z*1/N);
    for i=1:N
        y = y.*y1.*y2;
    end
end
function y = F2mf(z, beta, GM)
	y= sin((z+2).*atan((beta-1)./GM)).*(z+1);
    y1= (GM.^2+(beta-1).^2).^(-(z+2)/2/N);
	%y2= (GM).^(z*1/N);
    for i=1:N
        y = y.*y1.*y2;
    end
end
function y = F3pf(z, beta, GM)
	y= cos((z+3).*atan((beta+1)./GM)).*(z+1).*(z+2);
    y1= (GM.^2+(beta+1).^2).^(-(z+3)/2/N);
	%y2= (GM).^(z*1/N);
    for i=1:N
        y = y.*y1.*y2;
    end
end
function y = F3mf(z, beta, GM)
	y= cos((z+3).*atan((beta-1)./GM)).*(z+1).*(z+2);
    y1= (GM.^2+(beta-1).^2 ).^(-(z+3 )/2/N);
	%y2= (GM).^(z*1/N);
    for i=1:N
        y = y.*y1.*y2;
    end
end
function y = G1f(z, beta, GM)
	y= cos((z+1).*atan(2*beta./GM));
    y1= (GM.^2+4*beta.^2 ).^(-(z+1 )/2/N);
	%y2= (GM).^(z*1/N);
    for i=1:N
        y = y.*y1.*y2;
    end
end
function y = G2f(z, beta, GM)
	y= sin((z+2).*atan(2*beta./GM)).*(z+1);
    y1= (GM.^2+4*beta.^2 ).^(-(z+2 )/2/N);
	%y2= (GM).^(z*1/N);
    for i=1:N
        y = y.*y1.*y2;
    end
end
function y = G3f(z, beta, GM)
	y= cos((z+3).*atan(2*beta./GM)).*(z+1).*(z+2);
    y1= (GM.^2+4*beta.^2 ).^(-(z+3 )/2/N);
	%y2= (GM).^(z*1/N);
    for i=1:N
        y = y.*y1.*y2;
    end
end
function y = L1f(z, beta, GM)
	y= cos((z+1).*atan(2./GM));
    y1= (GM.^2+4 ).^(-(z+1 )/2/N);
	%y2= (GM).^(z*1/N);
    for i=1:N
        y = y.*y1.*y2;
    end
end
function y = L2f(z, beta, GM)
	y= sin((z+2).*atan(2./GM)).*(z+1);
    y1= (GM.^2+4 ).^(-(z+2 )/2/N);
	%y2= (GM).^(z*1/N);
    for i=1:N
        y = y.*y1.*y2;
    end
end
function y = L3f(z, beta, GM)
	y= cos((z+3).*atan(2./GM)).*(z+1).*(z+2);
    y1= (GM.^2+4 ).^(-(z+3)/2/N);
	%y2= (GM).^(z*1/N);
    for i=1:N
        y = y.*y1.*y2;
    end
end
function y = K1f(z, beta, GM)
    y=1;
	y1= (GM).^(-(z+1)/N);
	%y2= (GM).^(z*1/N);
    for i=1:N
        y = y.*y1.*y2;
    end
end
function y = K3f(z, beta, GM)
    y = (z+1).*(z+2);
	y1= (GM).^(-(z+3)/N);
	%y2= (GM).^(z*1/N);
    for i=1:N
        y = y.*y1.*y2;
    end
end
end