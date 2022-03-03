function Pq = saxs_starpolymer_benoit(q, Rg, f)
% H. Benoit, Journal of Polymer Science, 1953, 11, 507-510
v = Rg^2*f*q.^2/(3*f-2);
Pq = 2./(f*v.^2).*(v-1+exp(-v)+(f-1)/2.*(1-exp(-v)).^2);