function [F, V] = saxstruncatedoctahedron(qx, qy, qz, L)
% [F, V] = saxstruncatedoctahedron(qx, qy, qz, [L, b])
% L is the edge length of the original octahedron
% b is the edge length of the truncated pyramid
% truncated octahedron
% orientation is identical to the octahedron.
%
L0 = L*3; % edgelength of the original octahedron;
[F, V] = saxstruncatingoctahedron(qx, qy, qz, [L0, L]);