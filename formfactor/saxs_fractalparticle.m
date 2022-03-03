function y = saxs_fractalparticle(q, I0, avgR, Df)
% Sorensen & Wang, 1999, Phys. Rev. E Stat. Phys. Plasmas Fluids, 60, 7143.
y = sin((Df-1)*atan(q*avgR));
y = y./((Df-1).*q*avgR.*(1+(q*avgR).^2).^((Df-1)/2));
y = y*I0;