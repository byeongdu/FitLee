function Fcyl = saxsspiralcylinder(qx, qy, qz, parameter)
phi0 = parameter{1};
Ncyl = parameter{2};
cylgap = parameter{3};
Rcyl = parameter{4};
Lcyl = parameter{5};

j = sqrt(-1);
Fcyl = zeros(size(qx));
sizeq = size(qy);
for n=1:Ncyl;
    phi = phi0*(n-1);
    mat = [cos(phi*pi/180), -sin(phi*pi/180); sin(phi*pi/180), cos(phi*pi/180)];

    qyz = mat*[qy(:)';qz(:)'];
    qyt = qyz(1,:)';
    qzt = qyz(2,:)';
    qyt = reshape(qyt, sizeq);
    qzt = reshape(qzt, sizeq);
    F = saxscylinder(qx, qyt, qzt, [Rcyl, Lcyl]).*exp(-j*qx*(Rcyl*2+cylgap)*(n-1));
    Fcyl = Fcyl + F;
end