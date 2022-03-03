function F = saxslongRD(qx, qy, qz, param)
edgelength = param(1);
H = 1/sqrt(3)*edgelength;
R = sqrt(2)/sqrt(3)*edgelength;
L = 2/sqrt(3)*edgelength;
c = polyhedra('octant', R, L, H);
%cm = mean(c.vertices);
Forig = c.vertices(c.faces(2, :), :);
seq = planeconstruct(Forig, [0, 0, 0]);
Forig = Forig(seq, :);
a = polyhedra('longRD', param(1), param(2));
%figure(4);clf; drawpolyhedron(a);
Q = [qx(:), qy(:), qz(:)];
%RecT = [2*R, 2*R, 2H + param(2)];
RecT = [2*R, 2*R, param(2)];
rot = rotate_around_vector([0, 0, 1], 45);
Qn = Q*inv(rot);
Frec = saxscube(Qn(:,1),Qn(:,2),Qn(:,3),RecT/2);
F = Frec;
for i=1:numel(a.faces)
    if numel(a.faces{i}) == 4
        pf = a.vertices(a.faces{i}, :);
        t = abs(pf(:,3)) == param(1)/2;
        tp = (pf(:,1) == 0) & (pf(:,2)==0);
        A = pf(t, :);
        pf(t | tp, :) = [];
        tip = A; tip(3) = pf(1, 3);
        Ffin = [pf; A];
        seq = planeconstruct(Ffin, tip);
        Ffin = Ffin(seq, :);
        Ffin = Ffin - repmat(tip, 3, 1);
        invR = inv(Forig)*Ffin;
%        b = c;
%        b.vertices = b.vertices*invR+ repmat(tip, numel(b.vertices(:,1)), 1);
%        drawpolyhedron(b)
        Qn = Q*inv(invR);
        Fo = saxsoctant(Qn(:,1),Qn(:,2),Qn(:,3),[R, L, H]);
        F = F - Fo.*exp(sqrt(-1)*(tip(1)*qx+tip(2)*qy+tip(3)*qz));
    end
end
%c = polyhedra('pyramid', L, H);
%c.vertices = c.vertices + repmat([0, 0, 1]*H+[0, 0, param(2)/2], numel(c.vertices(:,1)), 1);
%c.color = 'r';
%drawpolyhedron(c)
Fpr = saxspyramid(qx, qy, qz, [L/2, H]);
F = F + Fpr.*exp(sqrt(-1)*qz*(H+param(2)/2));
Fpr = saxspyramid(qx, qy, -qz, [L/2, H]);
F = F + Fpr.*exp(-sqrt(-1)*qz*(H+param(2)/2));
