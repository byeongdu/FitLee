function Iq = saxs_cylinderavg(qxy, qz, funcname, parameters)
N = 10;
dtheta = pi/(2*N);
th = 0:dtheta:pi;
Iq = zeros(numel(qz), numel(qxy));
for k=1:numel(qz)
    [r, theta] = ndgrid(qxy, th);
    qx = r.*cos(theta);
    qy = r.*sin(theta);
    qzN = qz(k)*ones(size(qx));
    sizeqx = size(qx);

    qx = qx(:);
    qy = qy(:);
    qzN = qzN(:);
    tic
    F1 = feval(funcname, qx, qy, qzN, parameters);
    toc
    F1 = reshape(F1, sizeqx);
    y = abs(F1).^2;
    y = trapz(y, 2)*dtheta;
    Iq(k, :) = y';
end