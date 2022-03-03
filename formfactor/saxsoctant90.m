function F = saxsoctant90(qx,qy,qz,L)
i = sqrt(-1);
q1 = -qx*L/sqrt(2)+eps;
q2 = -qy*L/sqrt(2)+eps*2;
q3 = -qz*L/sqrt(2)+eps*3;
F = i*L^3/(2/sqrt(2));
F = F*(exp(i*q1)./(q1.*(q1-q2).*(q1-q3)+eps) + ...
    exp(i*q2)./(q2.*(q2-q1).*(q2-q3)+eps) + ...
    exp(i*q3)./(q3.*(q3-q1).*(q3-q2)+eps) - ...
    1./(q1.*q2.*q3+eps));