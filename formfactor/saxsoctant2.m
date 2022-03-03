function F = saxsoctant2(qx, qy, qz, RLH)
R = RLH(1);
L = RLH(2);
H = RLH(3);

tol = 1E-15;

qx = qx+eps;
qy = qy+eps;
qz = qz+eps;

%i = sqrt(-1);

%tanth = L/sqrt((2*R)^2-L^2);
%costh = sqrt(R^2-(L/2)^2)/R;
sinth = L/2/R;
th = asin(sinth);
costh = cos(th);
tanth = tan(th);

j = sqrt(-1);

%z = linspace(0, H, 20);
za = linspace(0, H, 50);
za = diff(za)/2+za(1:end-1);
dz = za(2)-za(1);
F = zeros(size(qx));
for i=1:numel(za)
    z = za(i);
    Rz = R - R/H*z;
    F0 = saxstriangle3(qx,qy,Rz, th);
    F = F + F0.*exp(-j*qz*z)*dz.*sinc(qz*dz/2);
end

%F = reshape(F, qxsize);
end
