function [F, V] = saxstruncatedoctahedron(qx, qy, qz, L, varargin)
% [F, V] = saxstruncatedoctahedron(qx, qy, qz, [L, b])
% L is the edge length of the original octahedron
% b is the edge length of the truncated pyramid
% truncated octahedron
% orientation is identical to the octahedron.
%
if numel(L) ==1
    L0 = L*3; % edgelength of the original octahedron;
else
    L0 = L(1);
    L = L(2);
end
[F, V] = saxstruncatingoctahedron(qx, qy, qz, [L0, L]);
if numel(varargin) == 2
    if strcmp(varargin{1}, 'shell thickness') 
        L02 = L0-varargin{2};
        L2 = L02/L0*L;
        [F2, V2] = saxstruncatingoctahedron(qx, qy, qz, [L02, L2]);
    end
    F = F - F2;
    V = V - V2;
end