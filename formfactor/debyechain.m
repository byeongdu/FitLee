function R = debyechain(vari, q)
% debyechain function
% ref : RJ, Roe book.
% eta = vari(1);
% v = vari(2);
% Rg = vari(3);

eta = vari(1);
v = vari(2);
Rg = vari(3);

x = Rg^2*q.^2;
D = 2*(exp(-1*x)+x-1)./x.^2;
R = eta^2*v^2*D;