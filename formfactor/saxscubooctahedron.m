function [F, V] = saxscubooctahedron(qx, qy, qz, L)
% y = saxscubooctahedron(qx, qy, qz, L)
% Edgelength : L
L0 = L*sqrt(2);
[F, V] = saxstruncatingcube(qx, qy, qz, [L0, L]);
