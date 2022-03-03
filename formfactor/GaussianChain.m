function test = GaussianChain(q,varia)
% this is for Gaussian Chain : GaussianChain(q, variable)
% variable is composed of 2 parameter(eta, v, n, l)
% cf n*l^2 = <Rg^2>
%  
%
eta = varia(1);
v = varia(2);
n = varia(3);
l = varia(4);
x = ((n*l^2)/6)*q.^2;
Dx = 2*(exp(-x)+x-1)./(x.^2);
test = eta^2*v^2*Dx;