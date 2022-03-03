function Z = paracrystal(q, a1, a2, a3, gfactor, N)
% function LatticeCal(qx, qy, qz, a1, a2, gfactor, N)
% this is the same function as BCCLattice.m ---- calculation for
% paracrystal
% but, it is generallized to every lattice, where a1, a2, a3 should be
% fundamental vectors of any lattice.
% this function calculate the paracrystal lattice factor with distortion of
% the second kind  ....... ref : Hashimoto et al. Macromolecules, Vol 27, 11,
% 1994, 3063
%
% finite size of domain is assumed in direction of a1 ==> N1
%                                                  a2 ==> N2 
%                                                  a3 ==> N3
% a1 and a2 should be a vector of cartesian coordinate...
% gfactor = [g11, g12, g21, g22];
% q could be matrix.
j = sqrt(-1);
qx = q(:,1);
qy = q(:,2);
qz = q(:,3);

if nargin<6
    N1 = -1;
    N2 = -1;
    N3 = -1;
    Ic1 = 0;
    Ic2 = 0;
    Ic3 = 0;
else
    N1 = N(1);
    N2 = N(2);
    N3 = N(3);
end

g = gfactor(1);

qa1 = a1(1)*qx + a1(2)*qy + a1(3)*qz;
qa2 = a2(1)*qx + a2(2)*qy + a2(3)*qz;
qa3 = a3(1)*qx + a3(2)*qy + a3(3)*qz;

%q = real(sqrt(qx.^2 + qy.^2 + qz.^2));
if sum(abs(a1)) > 0
    absF1 = abs(exp(-1/2*(g^2 * qa1.^2 + g^2 * qa2.^2 + g^2 * qa3.^2)));
    F1 = absF1.*exp(-j*qa1);
    Z1 = real((1 + F1)./(1 - F1));
    Z1(isnan(Z1)) = 0;
    Z1(isinf(Z1)) = 0;
    if N1>0
        Ic1 = -2*real(F1.*(1-F1.^N1)./((1-F1).^2));
    else
        Ic1 = 0;
    end
else
    Z1 = 1;
end

if sum(abs(a2)) > 0
    absF2 = abs(exp(-1/2*(g^2 * qa1.^2 + g^2 * qa2.^2 + g^2 * qa3.^2)));
    F2 = absF2.*exp(-j*qa2);
    Z2 = real((1 + F2)./(1 - F2));
    Z2(isnan(Z2)) = 0;
    Z2(isinf(Z2)) = 0;
    if N2>0
        Ic2 = -2*real(F2.*(1-F2.^N2)./((1-F2).^2 + eps));
    else
        Ic2 = 0;
    end
else
    Z2 = 1;
end

if sum(abs(a3)) > 0
    absF3 = abs(exp(-1/2*(g^2 * qa1.^2 + g^2 * qa2.^2 + g^2 * qa3.^2)));
    F3 = absF3.*exp(-j*qa3);
    Z3 = real((1 + F3)./(1 - F3));
    Z3(isnan(Z3)) = 0;
    Z3(isinf(Z3)) = 0;
    if N3>0
        Ic3 = -2*real(F3.*(1-F3.^N3)./((1-F3).^2 + eps));
    else
        Ic3 = 0;
    end
else
    Z3 = 1;
end

Z = (Z1 + Ic1/N1).*(Z2 + Ic2/N2).*(Z3 + Ic3/N3);

%Z = (Z1.*Z2.*Z3);