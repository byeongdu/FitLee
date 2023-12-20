function [Fq, Vtc] = saxstruncatedcube(qx, qy, qz, L, varargin)
% function F = saxstruncatedcube(qx, qy, qz, L)
% L = L0, where L0 is the edge length of a truncatedcube or
% L = [L0, f] for L0 is the length of a mother cube and f is degree of
% truncation.
% for an ideal TC, f = 0.4142 (or 1/(1+sqrt(2)))

if numel(L) == 1
    % truncatedcube, where fraction of truncating a cube is 0.4142. 
    L0 = L*(1+sqrt(2)); % edge length of the original cube.
else
    % truncating cube by a fraction (L(2))
    L0 = L(1);
    L = L0*L(2);
end
[Fq, Vtc]  = saxstruncatingcube(qx, qy, qz, [L0, L]);
if numel(varargin) == 2
    if strcmp(varargin{1}, 'shell thickness') 
        L02 = L0-varargin{2};
        L2 = L02/L0*L;
        [Fq2, Vtc2]  = saxstruncatingcube(qx, qy, qz, [L02, L2]);
    end
    Fq = Fq - Fq2;
    Vtc = Vtc - Vtc2;
end