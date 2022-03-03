function [F, V, Rmat, Rmatfaces, Rmatvertices, tri_tet] = saxsicosahedron( qx, qy, qz, L)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%a = 2;
V = 10*(3+sqrt(5))/3*(L/2)^3;

qx = qx + eps;
qy = qy + eps*2;
qz = qz + eps*3;


phi = (1+sqrt(5))/2;


vert = [0, 1, phi;
    0, -1, phi;
    0, 1, -phi;
    0, -1, -phi;
    1, phi, 0;
    -1, phi, 0;
    1, -phi, 0;
    -1, -phi, 0;
    phi, 0, 1;
    -phi, 0, 1;
    phi, 0, -1;
    -phi, 0, -1;
    ];
faces = [     1     2    10
     1    10     6
     1     6     5
     1     5     9
     1     9     2
     4     3    11
     4    11     7
     4     7     8
     4     8    12
     4    12     3
    12    10     8
    10     8     2
     8     2     7
     2     7     9
     7     9    11
     9    11     5
    11     5     3
     5     3     6
     3     6    12
     6    12    10];
t.vertices = vert;
t.faces = faces;

%tri_tet = [
%         1  -sqrt(3)  0
%      -1 -sqrt(3) 0
%      0 2/sqrt(3)  2*sqrt(2)/sqrt(3)
%];

v1 = t.vertices(t.faces(1, :), :);
edgeL = norm(v1(1,:)-v1(2,:));
Ru = distancePoints(t.vertices, [0,0,0]);% should be edgeL*sin(2*pi/5)
sizeqx = size(qx);
Qall = [qx(:), qy(:), qz(:)]';
F = zeros(size(Qall(1,:)));

%ico_tri = t.vertices(t.faces(1,:),:,:);
%plane = [ico_tri(1,:), ico_tri(1,:)-ico_tri(2,:), ico_tri(1,:)-ico_tri(3,:)];
%Htetra = distancePointPlane([0,0,0], plane)

% vertices of face 1
% [0 1 phi; 0 -1 phi; -phi 0 1];
% CM of face1 = [-phi/3, 0, (2*phi+1)/3];
% distance between [0,0,0] to CM of the face 1 = Htetra;
Htetra = 1/3*sqrt(phi^2 + (2*phi+1)^2);
H = Htetra/edgeL*L(1);


% when edgeL == 2
tri_tet = [
        1  -1/sqrt(3)  -Htetra
      -1 -1/sqrt(3) -Htetra
      0 2/sqrt(3)  -Htetra
];

inv_tri_tet = inv(tri_tet);




Rmat = {};
Rmatinv = {};
Rmatfaces = {};
Rmatvertices = {};
for i=1:length(t.faces)
    ico_tri = t.vertices(t.faces(i,:),:,:);
    %invM = tri_tet'*inv(ico_tri'); % for coordinate transform
    M = (inv_tri_tet*ico_tri)';
    invM = inv(M);
    
    %M = ico_tri\tri_tet;
    %M = M/abs(det(M))^(1/3);
    %invM = tri_tet\ico_tri;
    %invM = invM/abs(det(invM))^(1/3);
    
    if isempty(Rmat)
        Rmat{1} = M;
        Rmatinv{1} = invM;
        Rmatvertices{1} = (t.vertices(t.faces(1,:), :));
        Rmatfaces{1} = t.faces(1, :);
    else
        newrot = 1;
        for mati = 1:numel(Rmat);
            if isequal(sort(Rmatvertices{mati}), sort(-t.vertices(t.faces(i,:), :)))
                newrot = 0;
            end
        end
        if newrot
            n = numel(Rmat)+1;
            Rmat{n} = M;
            Rmatinv{n} = invM;
            Rmatvertices{n} = (t.vertices(t.faces(i,:), :));
            Rmatfaces{n} = t.faces(i, :);
        end
    end
end

for i=1:numel(Rmatinv);
    invM = Rmatinv{i};
    Q = invM*Qall;
    F = F + 2*real(saxsTet(Q(1,:), Q(2,:), Q(3,:), [L, H]));
%     M =Rmat{i};   % for particle transform.
%     temp = tri_tet*M';
%     plot3(temp(:,1), temp(:,2), temp(:,3), 'ro')
%     patch(temp(:,1), temp(:,2), temp(:,3), 'r')
end
% patch(tri_tet(:,1), tri_tet(:,2), tri_tet(:,3), 'g')
F = reshape(F, sizeqx);

function F = saxsTet(qx, qy, qz, LH)
    j = sqrt(-1);
    R = LH(1)/2;
    F = saxstetrahedron(qx, qy, qz, LH);
    F = F.*exp(j*qz*LH(2));% shift to -qz direction by LH(2)
end
end
