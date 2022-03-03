function [Fq, Vtc] = saxstruncatedcube(qx, qy, qz, L)
% function F = saxstruncatedcube(qx, qy, qz, L)
% L is the edge length of a cube.
L0 = L*(1+sqrt(2)); % edge length of the original cube.

[Fq, Vtc]  = saxstruncatingcube(qx, qy, qz, [L0, L]);