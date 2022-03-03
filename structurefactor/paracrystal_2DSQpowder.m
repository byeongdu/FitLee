function Za = paracrystal_2DSQpowder(qx, a1, gfactor, N)
% gfactor = da1/a1
% 2D Hexagonal paracrystal
qx = qx(:);
q = [qx, zeros(size(qx)), zeros(size(qx))];
ax = [a1, 0, 0];
ay = [0, a1, 0];
az = [0, 0, 0];
Za =  zeros(size(qx));
Nrot = 91;
ang = linspace(0, 90, Nrot);
for i=1:Nrot-1;
    R = rotate_around_vector([0, 0, 1], ang(i));
    R = R';
    a1 = ax*R;
    a2 = ay*R;
    a3 = az*R;
    if nargin == 4
        Z = paracrystal(q, a1, a2, a3, gfactor, N);
    else
        Z = paracrystal(q, a1, a2, a3, gfactor);
    end
    Za = Za + Z/(Nrot-1);
end