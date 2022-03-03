function p = logndistCDF2(r0, sig, x)
% p = schultzdistCDF(r0, sig, x);
% cumulative schultz distribution function
%p = normcdf(x, r0, sig);
p = zeros(size(x));
for i=1:numel(x)
%    dx = r0/20;
    dx = x(i)/100;
    t = 0:dx:x(i);
    k = logndist(r0, sig,t);
    p(i) = sum(k)*dx;
end
