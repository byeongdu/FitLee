function [y, name, pnames, pin] = ZhupairC(x, p, flag)
% [y, name, pnames, pin] = ZhuStruc(x, p, flag)
% Zhu et al. J. Chem. Phys., 1995, 104(22), 1996
% pnames=str2mat('a0','b', 'Df', 'Contrast^2', 'bkg');
if nargin == 2;
    w  = p(1);
    D = p(2);
%    r0 = p(3);
    %r = eps:1:50; r=r/50*0.5*D;lenr = length(r);
%    k = find(r<r0);r(k) = eps;
% ================== 1-d paracrystal... extended.
%Paracrystal = exp(pi*x.^2*w^2).*exp(-j*x.*D);
%In = (real(2*(1-Paracrystal.^50)./(1-Paracrystal)-1)).^2;y=In;
%=================================

%    y = 1-exp(-x/w).*sin(2*pi*x/D)./(2*pi*x/D);           %Zhu pair
%    correlation function

% ================== hard sphere,
%function SF = SFactor(q, Rhs, vf)
vf = w; Rhs = D;q = x;
alpha = (1 + 2*vf)^2/(1-vf)^4;
beta = -6*vf*(1 + vf/2)^2/(1-vf)^4;
gamma = vf*alpha/2;

if length(Rhs) > 1
    [R, q] = meshgrid(Rhs, q);
    A = R.*q*2;
else
    A = Rhs.*q*2;  
end
A = A+eps;
GA = alpha*(sin(A)-A.*cos(A))./A.^2 + ...
    beta*(2.*A.*sin(A) + (2 - A.^2).*cos(A)- 2)./A.^3 + ...
    gamma*(-1*A.^4.*cos(A) + 4*((3*A.^2-6).*cos(A) + (A.^3 - 6 * A).*sin(A) + 6))./A.^5;
y = 1./(1+24*vf*GA./A);
% ================================

else
	y=[];
	name='Zhu Structure factor';
	pnames=str2mat('w','D');
	if flag==1, pin=[10 100]; else pin = p; end
end