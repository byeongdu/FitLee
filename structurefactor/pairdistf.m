function Pdf_r = pairdistf(q_data, Intensity, real_coord)
% OK = pairdistf(q_data, Intensity, real_coord)
%

q = q_data;   

[xRc, yRc] = size(real_coord);
if xRc < yRc
    real_coord = real_coord';    % column vector
end

[xq, yq] = size(q_data);
if xq > yq
    q = q';                             % row vector.
end

[xI, yI] = size(Intensity);

if xI > yI
    Intensity = Intensity';         % row vector
end

temp = (real_coord*(Intensity.*q)).*sin(real_coord*q);
Pdf_r = sum(temp, 2)/2/pi.^2;