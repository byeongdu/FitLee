function [y, Rg] = saxs_poly_fractalparticle(q, I0, avgR, sig, Df)
% when primary nuclei are infintely small,
% the pair correlation function for an aggregate is
% gamma(r, Df, xsi) = r^(Df-3)h(r, xsi), 
% where h(r, xsi) are the cut-off function.
% When the cut off function is a Heaviside step function,
% h(r, xsi) = H(xsi - r) and polydispersity is introduced for xsi,
% where polydispersity function, w(xsi, avgR, z), is the schultz function.
% R. Besselink and J. E. ten Elshof, J. Appl. Cryst., 2014, 47, 1606.
z = (avgR/sig).^2 - 1;
Rg = (avgR/(z+1))*(1/2*Df*(Df+z+1)*(Df+z+2)/(Df+2))^(1/2);
if (Df <= 1) || (Df>3)
    error('Df should be equal or greater than 1 or less than 3')
end
zeta = floor(z);
SI = zeros(size(q));
for i=0:zeta
%   SI = SI + gamma(Df+i-1)/gamma(i+1)*sin((Df+i-1).*atan(q*avgR/(z+1)))...
%       ./(1+(q*avgR/(z+1)).^2).^((Df+i-1)/2);
    term1 = gamma(Df+i-1)./gamma(i+1);
    term2 = sin((Df+i-1).*atan(q*avgR/(z+1)))./(1+(q*avgR/(z+1)).^2).^((Df+i-1)/2);
    if isinf(term1) | isnan(term1)
        term1 = 1;
    end
    SI = SI + term1.*term2;
end

term1 = gamma(Df+zeta)/gamma(zeta+2);
if isinf(term1) | isnan(term1)
    term1 = 1;
end

SF = term1*sin((Df+zeta).*atan((q*avgR/(z+1))))...
    ./(1+(q*avgR/(z+1)).^2).^((Df+zeta)/2);
phi = z-zeta;

term1 = gamma(z+2)./gamma(Df+z+1);
if isinf(term1) | isnan(term1)
    term1 = 1;
end
y = I0*Df*term1./(q*avgR).*(SI+phi*SF);
