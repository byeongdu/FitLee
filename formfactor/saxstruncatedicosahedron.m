function saxstruncatedicosahedron(qx, qy, qz, L)

vert = [0, 1, 3*phi; 
    0, -1, 3*phi
    0, 1, -3*phi
    0, -1, -3*phi
    2, (1+2*phi), phi
    -2, (1+2*phi), phi
    2, -(1+2*phi), phi
    -2, -(1+2*phi), phi
    2, (1+2*phi), -phi
    -2, (1+2*phi), -phi
    2, -(1+2*phi), -phi
    -2, -(1+2*phi), -phi
    1, (2+phi), 2*phi
    -1, (2+phi), 2*phi
    1, -(2+phi), 2*phi
    -1, -(2+phi), 2*phi
    1, (2+phi), -2*phi
    -1, (2+phi), -2*phi
    1, -(2+phi), -2*phi
    -1, -(2+phi), -2*phi
];
vert=[vert;[vert(:,2),vert(:,3), vert(:,1)];[vert(:,3), vert(:,1), vert(:,2)]];
F = findfaces(vert);
