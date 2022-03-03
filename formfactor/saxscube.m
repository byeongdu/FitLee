function [F, V] = saxscube(qx, qy, qz, parameter)
% Form factor of Cube or rectangular parallelpipedons
% F = saxscube(qx, qy, qz, parameter)
% parameter is a or [a, b, c]
% variable, x, p, flag
% cube: length of each cell is 2*a (or 2*b, 2*c)...

if numel(parameter) == 1
    a = parameter(1);
    b = parameter(1);
    c = parameter(1);
    V = 8*a^3;
else
    a = parameter(:,1)';
    [~,ny] = size(parameter);
    if (ny == 1) 
        b = a; c = a;
    elseif (ny == 2)
        b = parameter(:,2)';
        c = a;
    elseif (ny == 3)
        b = parameter(:,2)';
        c = parameter(:,3)';
    end
    V = 8.*a.*b.*c;
end
F = sinc(qx*a).*sinc(qy*b).*sinc(qz*c)*V;