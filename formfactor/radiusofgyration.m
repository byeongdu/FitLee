function Rg = radiusofgyration(p)
% p is a coordinates of atoms (=[x0,y0,z0; x1, y1, z1;...;xn,yn,zn]

if numel(p)<3000
    [p0x, p1x] = meshgrid(p(:,1), p(:,1));
    [p0y, p1y] = meshgrid(p(:,2), p(:,2));
    [p0z, p1z] = meshgrid(p(:,3), p(:,3));

    dsqr = (p0x-p1x).^2+(p0y-p1y).^2+(p0z-p1z).^2;
    Rgsqr = sum(sum(dsqr))/2/numel(dsqr);
else
    [nx, ny] = size(p);
    dsqr = 0;
    for i=1:nx
        k = p(:,:)-repmat(p(i,:), nx, 1);
        dsqr = dsqr + sum(sum(k(:,1).^2+k(:,2).^2+k(:,3).^2));
    end
    Rgsqr = dsqr/2/nx^2;
end
Rg = sqrt(Rgsqr);