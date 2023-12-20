function [F, V] = saxscubooctahedron(qx, qy, qz, L, varargin)
% y = saxscubooctahedron(qx, qy, qz, L)
% Edgelength : L
L0 = L*sqrt(2);

[F, V] = saxstruncatingcube(qx, qy, qz, [L0, L]);

if numel(varargin)==2
    if strcmp(varargin{1}, 'shell thickness')
        sht = varargin{2};
    end
    L = L-sht;
    L0 = L*sqrt(2);
    [F2, V2] = saxstruncatingcube(qx, qy, qz, [L0, L]);
    F = F-F2;
    V = V-V2;
end
