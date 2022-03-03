function F = saxspolygonalcon(qx, qy, qz, NLH)
N = NLH(1);
L = NLH(2);
H = NLH(3);
Qall = [qx(:), qy(:), qz(:)];
F = zeros(size(qx(:)));
R = L/2/sin(pi/N);
for i=1:N
    a = 2*pi/N*i;
    Mt = [cos(a), sin(a), 0;...
        -sin(a), cos(a), 0;...
        0, 0, 1];
    Q = Qall*Mt';
    Ft = saxsoctant(Q(:,1), Q(:,2), Q(:,3), [R, L, H]);
    F = F+Ft;
end