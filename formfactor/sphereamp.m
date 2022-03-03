function [Formfactor, q, R] = sphereamp(q, r, sig, flag)
% this Amplitude of Hard sphere type particle : spheretype(q, r, flaq)
% this function is valid in the region 
% When r is a scaler, formfactor matrix will be an array with size of numel(q)
% when r is a array,  formfactor matrix will be a 2D matrix with size of numel(q)
% q should be raw vector
% p = [R, sig, I0]
% sig is interface roughness
dims = size(r);
if prod(dims) == max(dims)
    dims = max(dims);
end
%q = q(:);
numq = numel(q);
k = find(abs(q) < eps);
if ~isempty(k)
%    if (k>0) & (k<length(q))
        q(k) = 1E-8;
        %   end
end
V = 4*pi/3*r.^3;
if nargin == 2
    if size(r) == 1
        %test = eta^2*V^2*9*((sin(R.*q) - R.*q.*(cos(R.*q))).^2)./((R.*q).^6).*exp(-1.*sig.^2*q.^2);
        R = r;
        test = 3*V*(sin(R.*q) - R.*q.*(cos(R.*q)))./((R.*q+eps).^3);
    else
        [R, q] = meshgrid(r, q);
        R = reshape(R,[numq, dims]);
        q = reshape(q,[numq, dims]);
        %[R, q] = ndgrid(r, q);
        %R = reshape(R,[dims, numq]);
        %q = reshape(q,[dims, numq]);

        V = 4*pi/3*R.^3;
        test = 3*V.*(sin(R.*q) - R.*q.*(cos(R.*q)))./((R.*q+eps).^3);
        %test = eta.^2.*V.^2.*9.*((sin(R.*q) - R.*q.*(cos(R.*q))).^2)./((R.*q).^6).*exp(-1.*sig.^2.*q.^2);
    end
    Formfactor = test;
    %Formfactor = sum(test, 2);
    %Formfactor = reshape(test, dims);
    %Formfactor = test.*exp(-1/2*sig^2*q.^2);
elseif nargin == 3
    if size(r) == 1
        %test = eta^2*V^2*9*((sin(R.*q) - R.*q.*(cos(R.*q))).^2)./((R.*q).^6).*exp(-1.*sig.^2*q.^2);
        R = r;
        test = 3*V*(sin(R.*q) - R.*q.*(cos(R.*q)))./((R.*q+eps).^3);
    else
        [R, q] = meshgrid(r, q);
        test = 3*V*(sin(R.*q) - R.*q.*(cos(R.*q)))./((R.*q+eps).^3);
        %test = eta.^2.*V.^2.*9.*((sin(R.*q) - R.*q.*(cos(R.*q))).^2)./((R.*q).^6).*exp(-1.*sig.^2.*q.^2);
    end
    Formfactor = sum(test, 2);
    %Formfactor = test.*exp(-1/2*sig^2*q.^2);
    %Formfactor = reshape(Formfactor, dims);
else
   Formfactor=[];
   name='HardSphere Amplitude';
   
    pnames=str2mat('R');
	if flag==1, pin=[10]; else pin = p; end
end