function [F, V] = saxsstellatetetrahedron(qx, qy, qz, L)
% [F, V] = saxsstellatetetrahedron( qx, qy, qz, L)
% L should be [edge length of icosahedron, Height for stellate]
% if you like a regular tetrahedron, h = sqrt(2/3)*L
%
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%Hstellate = sqrt(2/3)*L;
Hstellate = L(1);
if numel(L) == 1
    Hstellate = sqrt(2/3)*L;
end
if numel(L) == 2
    %h = L(2);
    Hstellate = L(2); %Hstellate*h;
end
L = L(1);

%qx = qx + eps;
%qy = qy + eps;
%qz = qz + eps;

%a = 2;
t = polyhedra('tetrahedron', 2);
cmz = sqrt(6)/12;% center of mass of tetrahedron is at sqrt(6)/12*L above the bottom plane center.
cm = [0,0,cmz*2]; 
j = sqrt(-1);
Ftet = saxstetrahedron(qx,qy,qz,L(1));
%t2 = t;t2.vertices = t2.vertices/2*L;
%Ftet = saxsPolyhedronAmp(qx, qy, qz, t2);
Ftet = Ftet.*exp(-j*qz*cmz*(-L));



sizeqx = size(qx);
Qall = [qx(:), qy(:), qz(:)]';

%SHtet = Hstellate/L*2;

t.vertices = t.vertices-repmat(cm, [4,1]);
%bodytet = t.vertices(t.faces(1,:),:,:);
%plane = [bodytet(1,:), bodytet(1,:)-bodytet(2,:), bodytet(1,:)-bodytet(3,:)];
%Htetra = distancePointPlane([0,0,0], plane);

stellate_tri_tet = [
         1  -1/sqrt(3) -cmz*2
      -1 -1/sqrt(3) -cmz*2
      0 2/sqrt(3)  -cmz*2
];

% or, you can do..
% stellate_tri_tet = t.vertices(t.faces(1,:),:,:);

%tstellate = polyhedra('tetrahedron', 2, SHtet);
%tstellate.vertices = tstellate.vertices - repmat([0,0,SHtet], 4, 1);
V = t.vertices;
F = t.faces;
NV = size(F, 2);
NF = length(F);
%figure;hold on;
Fs = zeros(size(Qall(1,:)));

for i=1:NF;
    bodytet = V(F(i,:),:,:);
    % Now let's start to calculate F of stellate tetrahedron
    %plane = [bodytet(1,:), bodytet(1,:)-bodytet(2,:), bodytet(1,:)-bodytet(3,:)];
    %D = distancePointPlane([0,0,0], plane);
    %n(i,:) = n(i,:)/norm(n(i,:));
    %D = abs(D);
    %mirror_ico_tri = bodytet - (D+SHtet)*repmat(n(i,:), [3,1]);
    %invM_for_stellate = stellate_tri_tet'*inv(mirror_ico_tri');
    invM_for_stellate = stellate_tri_tet'*inv(bodytet');
    Q = invM_for_stellate*Qall;

%    M = inv(invM_for_stellate);
%    tn = (tstellate.vertices)*M'+(D+SHtet)*repmat(n(i,:), [4,1]);
%    drawpolyhedron(tn*L/2, tstellate.faces, rand(3,1));
    Fs = Fs + saxsStellate(Q(1,:), Q(2,:), -Q(3,:), [L, Hstellate], cmz*L);
end
%drawpolyhedron(t.vertices*L/2, t.faces, rand(3,1));
F = Ftet + reshape(Fs, sizeqx);
V=0;
end

function F = saxsStellate(qx, qy, qz, LH, Hshift)
    % Hshift should be distance between [0,0,0] and a plane of
    % icosa_triangle..
    j = sqrt(-1);
%    R = LH(1)/2;
    %tstellate = polyhedra('tetrahedron', LH(1), LH(2));
    %F = saxsPolyhedronAmp(qx, qy, qz, tstellate);
    F = saxstetrahedron(qx, qy, qz, LH);
    F = F.*exp(-j*qz*Hshift);% shift to qz direction by LH(2)
end
