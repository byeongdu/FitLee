function [F, V] = saxsstellateicosahedron( qx, qy, qz, L)
% [F, V] = saxsstellateicosahedron( qx, qy, qz, L)
% L should be [edge length of icosahedron, h for stellate]
% when h=1, height of tetrahedron of the stallate is 1
% if you like a regular tetrahedron, h = sqrt(2/3)
%
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Hstellate = L(1);
if numel(L) == 1
    Hstellate = sqrt(2/3)*L;
end
if numel(L) == 2
    %h = L(2);
    Hstellate = L(2); %Hstellate*h;
end
L = L(1);

%a = 2;
V = 10*(3+sqrt(5))/3*(L/2)^3;
[~,Vt]=saxstetrahedron(0,0,0,[L,Hstellate]);
V = V+Vt*20;% 20 stellate tetrahedron..

qx = qx + eps;
qy = qy + eps;
qz = qz + eps;


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

%edgeL = norm(t.vertices(1,:)-t.vertices(2,:)); % edge length should be 2
%Ru = distancePoints(t.vertices, [0,0,0]);% should be edgeL*sin(2*pi/5)
sizeqx = size(qx);
Qall = [qx(:), qy(:), qz(:)]';
F = zeros(size(Qall(1,:)));
Fs = zeros(size(Qall(1,:)));

ico_tri = t.vertices(t.faces(1,:),:,:);
plane = [ico_tri(1,:), ico_tri(1,:)-ico_tri(2,:), ico_tri(1,:)-ico_tri(3,:)];
Htetra = distancePointPlane([0,0,0], plane);
H = Htetra/2*L(1);

tri_tet = [
         1  -1/sqrt(3)  -Htetra
      -1 -1/sqrt(3) -Htetra
      0 2/sqrt(3)  -Htetra
];

SHtet = Hstellate/L*2;
stellate_tri_tet = [
         1  -1/sqrt(3)  -SHtet
      -1 -1/sqrt(3) -SHtet
      0 2/sqrt(3)  -SHtet
];


%figure;hold on;
for i=1:length(t.faces)
    ico_tri = t.vertices(t.faces(i,:),:,:);
    invM = tri_tet'*inv(ico_tri'); % for coordinate transform
    Q = invM*Qall;
    F = F+saxsTet(Q(1,:), Q(2,:), Q(3,:), [L, H]);
    % Now let's start to calculate F of stellate tetrahedron
    % 
    %cm = mean(ico_tri);
    %n = -cm;
    %mirror_ico_tri = ico_tri + (1+Hstellate/H)*repmat(n, [3,1]);
    %invM_for_stellate = stellate_tri_tet'*inv(mirror_ico_tri');
    
    %Q = invM_for_stellate*Qall;
    
    % Don't need to reflect the plane.
    % Stellate is on the opposite side of the tetrahedron making icosa.
    % therefore the same rotation matrix can be used if only, sign of qz is
    % changed later on...
    Fs = Fs+saxsStellate(Q(1,:), Q(2,:), -Q(3,:), [L, Hstellate], H);
end
F = F + Fs;
F = reshape(F, sizeqx);

function F = saxsTet(qx, qy, qz, LH)
    j = sqrt(-1);
    %R = LH(1)/2;
    F = saxstetrahedron(qx, qy, qz, LH);
    F = F.*exp(j*qz*LH(2));% shift to -qz direction by LH(2)
end
function F = saxsStellate(qx, qy, qz, LH, Hshift)
    % Hshift should be distance between [0,0,0] and a plane of
    % icosa_triangle..
    j = sqrt(-1);
    %R = LH(1)/2;
    F = saxstetrahedron(qx, qy, qz, LH);
    F = F.*exp(-j*qz*Hshift);% shift to qz direction by LH(2)
end
end
