function [Formfactor, V2, V] = spherePq(q, r, flag)
% this is for Hard sphere type particle : spheretype(q, r, flaq)
% this function is valid in the region 
% q should be raw vector

if nargin == 2
    eta = 1;
    V = (4/3*pi*r.^3);
    V2 = V.^2;
    if size(r) == 1
        %test = eta^2*V^2*9*((sin(R.*q) - R.*q.*(cos(R.*q))).^2)./((R.*q).^6).*exp(-1.*sig.^2*q.^2);
        R = r;
        test = eta^2*V2*9*((sin(R.*q) - R.*q.*(cos(R.*q))).^2)./((R.*q).^6);
    else
        V = repmat(V, [numel(q), 1]);
        [R, q] = meshgrid(r, q);
        test = eta.^2.*V.^2.*9.*((sin(R.*q) - R.*q.*(cos(R.*q))).^2)./((R.*q).^6);
        %test = eta.^2.*V.^2.*9.*((sin(R.*q) - R.*q.*(cos(R.*q))).^2)./((R.*q).^6).*exp(-1.*sig.^2.*q.^2);
    end
    %Formfactor = sum(test, 2);
    Formfactor = test;
t = abs(q)< eps;
Formfactor(t) = V2;
else
    
   Formfactor=[];
   name='HardSphere Model';
   
    pnames=str2mat('R');
	if flag==1, pin=[10]; else pin = p; end
end       