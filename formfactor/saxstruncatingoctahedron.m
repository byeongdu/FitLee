function [F, V] = saxstruncatingoctahedron(qx, qy, qz, Lb)
% [F, V] = saxstruncatedoctahedron(qx, qy, qz, [L, b])
% L is the edge length of the original octahedron
% b is the edge length of the truncated pyramid
% truncated octahedron
% orientation is identical to the octahedron.
%
j = sqrt(-1);
%if numel(parameter) == 1
%    L = parameter(1)*2;
%    b = L/3;
%elseif numel(parameter) == 2
%    L = parameter(1)*2;
%    b = parameter(2);
%end

if size(Lb,2) == 1
    Lm = Lb(:,1);
    bm =Lm/3;   % equal edge length pyramid...
else
    Lm = Lb(:,1);bm = Lb(:,2);
end

Q = [qx(:), qy(:), qz(:)];

if numel(Lm > 1)
    [Ln, ~] = meshgrid(Lm, qx);
    [bn, ~] = meshgrid(bm, qy);
    [~, qzn] = meshgrid(Lm, qz);
    L = Lm;
    b = bm;
    sizeqzn = size(qzn);
else
    L = Lm;Ln=L;
    b = bm;bn=b;
    sizeqzn = size(qz);
end


[F, Voch] = saxsoctahedron(qx, qy, qz, L/2);
if b == 0
    V = Voch;
    return
end
Hto = (Ln-bn)/sqrt(2); % [numel(qx) by numel(L)] matrix or a scalar
[Fpy1, Vp] = saxspyramid(qx, qy, qz, b/2);
Fpy1 = Fpy1.*exp(-j*qzn.*Hto);
Fpy2 = saxspyramid(qx, qy, -qz, b/2);
Fpy2 = Fpy2.*exp(j*qzn.*Hto);
ang = [pi/2, -pi/2];

for i=1:2
    ct = cos(ang(i));
    st = sin(ang(i));
    R1 = [(1+ct)/2, (1-ct)/2, 1/sqrt(2)*st;
        (1-ct)/2, (1+ct)/2, -1/sqrt(2)*st;
        -1/sqrt(2)*st, 1/sqrt(2)*st, ct];
    R2 = [(1+ct)/2, -(1-ct)/2, -1/sqrt(2)*st;
        -(1-ct)/2, (1+ct)/2, -1/sqrt(2)*st;
        1/sqrt(2)*st, 1/sqrt(2)*st, ct];
    Qrot1 = Q*R1;
    Qrot2 = Q*R2;
    if numel(Hto)>1
        [~, qznrot1] = meshgrid(Lm, Qrot1(:,3));
        [~, qznrot2] = meshgrid(Lm, Qrot2(:,3));
    else
        qznrot1 = Qrot1(:,3);
        qznrot2 = Qrot2(:,3);
    end
    F1=saxspyramid(Qrot1(:,1),Qrot1(:,2),-Qrot1(:,3),b/2).*exp(j*qznrot1.*Hto);
    Fpy1 = Fpy1 + reshape(F1, sizeqzn);
    F2=saxspyramid(Qrot2(:,1),Qrot2(:,2),-Qrot2(:,3),b/2).*exp(j*qznrot2.*Hto);
    Fpy2 = Fpy2 + reshape(F2, sizeqzn);
end
F = F-Fpy1-Fpy2;
V = Voch-Vp*6;