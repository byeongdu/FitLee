function [Rmat, Rmatinv, Rmatfaces, Rmatvertices, Vbase, H] = getRotMat_Subunit2Face(vertices, facesofinterest)
% [Rmat, Rmatinv, Rmatfaces, Rmatvertices, Vbase, H] = getRotMat_Subunit2Face(vertices, facesofinterest)
% This function calculate R and inverse R mat for a polyhedron.
% R matrix is a matrix to rotate the subunit to a face of the polyhedron.
% INPUT:
%   vertices: coordinates of all vertices. N_vert x 3
%   faces_of_interest: faces (set of indices of vertices). N_face x N_gon
%       This should not be a cell. So if a soccerball is of interest, 
%       decide whether you want to calculate for pentagonal faces or
%       hexagonal faces.
% OUTPUT:
%   R
%   inverse_R
%   R_faces
%   R_vertices
%   Vbase : vertices of the base of the subunit.
%       you can check the orientation of the polyhedron.
%       Rotated_Vertices = Vbase * inverse_R;
%   H: if you multiply L to H then you will get the 

obj.vertices = vertices;
obj.faces = facesofinterest;
[Nfaces, N] = size(facesofinterest); % N-polygonal cone, N : number of vertices in a face.

v1 = obj.vertices(obj.faces(1, :), :);

edgeL = norm(v1(1,:)-v1(2,:));
ico_tri = obj.vertices(obj.faces(1,:),:,:);
plane = [ico_tri(1,:), ico_tri(1,:)-ico_tri(2,:), ico_tri(1,:)-ico_tri(3,:)];
Htetra = abs(distancePointPlane([0,0,0], plane));

subunitbase = [
        edgeL/2  -edgeL/2/tan(pi/N)  -Htetra
      -edgeL/2 -edgeL/2/tan(pi/N) -Htetra
      0 0  -Htetra
];

inv_tri_tet = inv(subunitbase);

Vbase = zeros(N, 3);
for i=1:N
    Rpolygon = [cos(2*pi/N*i), sin(2*pi/N*i), 0;...
        -sin(2*pi/N*i), cos(2*pi/N*i), 0;...
        0, 0, 1];
    Vbase_tmp = subunitbase(1, :)*Rpolygon;
    Vbase(i, :) = Vbase_tmp;
end
H = Htetra/edgeL;



Rmat = {};
Rmatinv = {};
Rmatfaces = {};
Rmatvertices = {};
for i=1:Nfaces
    tri = obj.vertices(obj.faces(i,:),:,:);
    CM = mean(tri);two_nearest_vertices=[];
    for k=2:numel(tri(:,1));
        if norm(tri(k,:)-tri(1,:)) <= edgeL+edgeL/10;
            two_nearest_vertices = [tri(1,:); tri(k,:)];
            break
        end
    end
    ico_tri = [two_nearest_vertices; CM];
    
    invM = inv_tri_tet*ico_tri;
    M = inv(invM);
    
    
    if isempty(Rmat)
        Rmat{1} = M;
        Rmatinv{1} = invM;
        Rmatvertices{1} = (obj.vertices(obj.faces(1,:), :));
        Rmatfaces{1} = obj.faces(1, :);
    else
        newrot = 1;
        for mati = 1:numel(Rmat);
            if isequal(sort(Rmatvertices{mati}), sort(-obj.vertices(obj.faces(i,:), :)))
                newrot = 0;
            end
        end
        if newrot
            n = numel(Rmat)+1;
            Rmat{n} = M;
            Rmatinv{n} = invM;
            Rmatvertices{n} = (obj.vertices(obj.faces(i,:), :));
            Rmatfaces{n} = obj.faces(i, :);
        end
    end
end