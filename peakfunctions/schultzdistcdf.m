function y = schultzdistcdf(r0, sig, x)
% p = schultzdistCDF(r0, sig, x);
% cumulative schultz distribution function
p = zeros(size(x));
%maxd = schultzdist(1E5, r0, sig)
%for i=1:numel(x)
%    dx = r0/20;
%    t = 0:dx:x(i);
%    k = schultzdist(t, r0, sig);
%    p(i) = sum(k)*dx;
%end
%dx = sig/r0;
%t = 0.0:dx:x(end);
if numel(x) == 1
    if x ==inf
        dx = sig/r0*10;
        if dx < 10
            t = 0:dx:r0*10;
        else
            t = 0:dx:1E5;
        end
        k = schultzdist(t, r0, sig);
        p = cumsum(k)*dx;
        y = p(end);
    return
    else
        dx = sig/r0*10;
        t=0:dx:x(end);
        k = schultzdist(t, r0, sig);
        p = cumsum(k)*dx;
        y = p(end);
    end
else
    dx = sig/r0*10;
    if dx>x(end)
        dx = x(end)/2;
    end
    t=0:dx:x(end);
    k = schultzdist(t, r0, sig);
    p = cumsum(k)*dx;
    y = interp1(t,p,x);
end
y = y/schultzdistcdf(r0, sig, inf);
%z = (r0/sig).^2-1;
%p = ((z+1)/r0)^(z+1)/z*r0/(z+1).*x.^z.*(z/r0.*x).^(-1/2*z).*exp(-1/2*z/r0.*x).*WhittakerM(1/2*z, 1/2+1/2*z, z/r0*x)/gamma(z+1);
%end
%function y = WhittakerM(k, m, z)
%y = z.^(1/2+m).*exp(-z/2).*(1+(1/2+m-k)/(1*(2*m+1)).*z + (1/2+m-k)*(3/2+m-k)/(2*(2*m+1)*(2*m+2))*z.^2);
%end