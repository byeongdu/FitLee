function y=SchultzCoreShellFF3(q, p)
% Analytical solution for a polydisperse spherical core-shell particles.
% Only core is polydisperse
% Numeric

    Rcore  = p(1);
    Rshell = p(2);
    Rcore_width = p(3);
    rho1 = p(4);
    rho2 = p(5);
    rho3 = p(6);

    beta = Rshell/Rcore;

    [gDist, xx] = logndist9(Rcore, Rcore_width);

    [x1, x2]= size(q);
    
    if x1 < x2
        q = q';
    end

    y = zeros(size(q));
    Norm = 0;
    for i = 1:length(xx)
	  Rcore = xx(i);
	  Rshell = beta*Rcore;
        Formfactor = (rho1-rho2)*(4/3*pi*Rcore^3)*sphereAmp1(q, Rcore) + (rho2-rho3)*(4/3*pi*Rshell^3)*sphereAmp1(q, Rshell);
	  Normfactor = (rho1-rho2)*(4/3*pi*Rcore^3) + (rho2-rho3)*(4/3*pi*Rshell^3);
        y = y + gDist(i)*abs(Formfactor).^2;
	  Norm = Norm + gDist(i)*Normfactor^2;
    end
    %y = y/Norm;

function test = sphereAmp1(q, r)
test = 3*(sin(q*r) - r*q.*(cos(q*r)))./((q*r).^3);
