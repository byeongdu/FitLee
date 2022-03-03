function [F, V] = saxspentagonalbipyramid(qx, qy, qz, LH)

[F1, V1] = saxspentagonalcone(qx, qy, qz, [LH(:,1), LH(:,2)/2]);
[F2, V2] = saxspentagonalcone(qx, qy, -qz, [LH(:, 1), LH(:,2)/2]);
F = F1+F2;
V = V1+V2;