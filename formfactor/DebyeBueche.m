function [y, name, pnames, pin]=DebyeBueche(x, p, flag)
% DebyeBueche(q, vari)
% I0 = vari(1);
% xi = vari(2);
% C = vari(3) : compressibility effect.

if nargin==2
    if numel(p)==3
        I0 = p(1);
        xi = p(2);
        C = p(3);
    elseif numel(p)==4
        % p(1) = delta_rho;
        % p(2) = volume fraction
        % p(3) = xsi
        % p(4) = constant background
        % This will return absolute intensity.
        xi = p(3);
        C = p(4);
        I0 = 8*pi*p(1)^2*p(2)*xi^3;
    end
    y = I0./(1+xi.^2.*x.^2).^2+C;

else
	y=[];
	name='DebyeBueche model';
	pnames=str2mat('I0', 'xi', 'Background');
	if flag==1 
        pin=[1 130 0]; 
    else 
        pin = p; 
    end
end
