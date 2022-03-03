function [Formfactor, name, pnames, pin] = spheretype(q, r, flag)
% this is for Hard sphere type particle : spheretype(q, r, flaq)
% this function is valid in the region 
% q should be raw vector

if nargin == 2
    eta = 1;
    V = 1;
    q(find(q==0))=eps;
    r(find(r==0))=eps;
    if size(r) == 1
        %test = eta^2*V^2*9*((sin(R.*q) - R.*q.*(cos(R.*q))).^2)./((R.*q).^6).*exp(-1.*sig.^2*q.^2);
        R = r;
        test = 9*((sin(R.*q) - R.*q.*(cos(R.*q))).^2 + eps)./((R.*q).^6 + eps);
    else
        [R, q] = meshgrid(r, q);
        test = eta.^2.*V.^2.*9.*((sin(R.*q) - R.*q.*(cos(R.*q))).^2)./((R.*q).^6);
        %test = eta.^2.*V.^2.*9.*((sin(R.*q) - R.*q.*(cos(R.*q))).^2)./((R.*q).^6).*exp(-1.*sig.^2.*q.^2);
    end
    %Formfactor = sum(test, 2);
    Formfactor = test;

else
    
   Formfactor=[];
   name='HardSphere Model';
   
    pnames=str2mat('R');
	if flag==1, pin=[10]; else pin = p; end
end       