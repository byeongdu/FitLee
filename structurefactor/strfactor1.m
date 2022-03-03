function [y, name, pnames, pin] = strfactor1(q, p, flag)
% structure factor under gaussian potential
% Gaussian core model
% Ref : Macromolecules, 2001, 34, 9, 2914.
%       PR, E, 2000, 62, 7961
%       J. Phys: Condens. Matter, 2000, 12, 5087
%    be = p(1);
%    ep = p(2);
%    rho = p(3);
%    sig = p(4);

if nargin==2;

    be = p(1);
    ep = p(2);
    rho = p(3);
    sig = p(4);
    y = 1./(1+ (pi).^(3/2)*be*ep*rho*sig.^3*exp(-(q*sig/2).^2));
else
	y=[];
	name='GaussianCoreModel_Structurefactor';
	pnames=str2mat('Beta','Epsion','Rho','Sigma');
	if flag==1, pin=[0.6 0.6 0.6 20]; else pin = p; end
end
